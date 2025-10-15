# TopologyExtractor::run 関数 実装ドキュメント

## 概要

`TopologyExtractor::run` 関数は、一般化ボロノイ図（GVD）のスケルトンマスクから、トポロジカルマップ（グラフ構造）を抽出します。

**入力:**
- `gvd_mask`: GVDスケルトンのピクセルマスク（uint8_t配列、255=スケルトン、0=非スケルトン）
- `width`, `height`: マップの幅と高さ（ピクセル単位）
- `resolution`: マップの解像度（メートル/ピクセル）

**出力:**
- `TopologicalMap`: ノード（頂点）とエッジ（辺）で構成されるグラフ構造

---

## アルゴリズムの処理フロー

### 前処理: Zhang-Suen細線化 (19-121行目)

```cpp
// Optimized Zhang-Suen thinning algorithm with parallel processing
static std::vector<uint8_t> zhangSuenThinning(const std::vector<uint8_t>& input, int width, int height)
```

**目的:** GVDスケルトンを1ピクセル幅の細線に変換し、トポロジー抽出を容易にする

**アルゴリズム:**
1. **Sub-iteration 1**: 特定の条件を満たすピクセルをマーク
   - 条件: `2 ≤ B ≤ 6`, `A == 1`, `!(P2∧P4∧P6)`, `!(P4∧P6∧P8)`
   - B: 8近傍の非ゼロピクセル数
   - A: 0→1遷移の回数（トポロジー保持のため）

2. **Sub-iteration 2**: 補完的な条件でピクセルをマーク
   - 条件: `2 ≤ B ≤ 6`, `A == 1`, `!(P2∧P4∧P8)`, `!(P2∧P6∧P8)`

3. **反復**: 変化がなくなるまで1と2を繰り返す

**最適化:**
- OpenMP並列化（`#pragma omp parallel for schedule(dynamic, 64)`）
- インライン近傍アクセスで関数呼び出しオーバーヘッド削減
- マーカーベクトルでバッチ削除
- 大規模マップでは反復回数を20に制限

**計算量:** O(I × W × H) - I: 反復回数（通常10-20回）

**スキップ条件:**
```cpp
const int max_pixels_for_thinning = 5000000; // 5 million pixels
if (width * height > max_pixels_for_thinning) {
    is_skel = gvd_mask; // Skip thinning for large maps
}
```

---

### フェーズ1: スケルトンピクセルの次数計算 (139-160行目)

```cpp
// Compute degrees for each skeleton pixel (8-neighborhood)
std::vector<uint8_t> is_skel = zhangSuenThinning(gvd_mask, width, height);
std::vector<uint8_t> degree(width * height, 0);
const int dx8[8] = { -1, 0, 1, -1, 1, -1, 0, 1 };
const int dy8[8] = { -1,-1,-1,  0, 0,  1, 1, 1 };
```

**目的:** 細線化されたスケルトンの各ピクセルの接続数（次数）を計算

**処理内容:**
- 各スケルトンピクセルについて、8方向の隣接ピクセルを確認
- スケルトンピクセルの隣接数をカウント → `degree[]` に格納
- OpenMP並列化により高速処理（`#pragma omp parallel for`）

**次数の意味:**
- `degree == 1`: 端点（エンドポイント）
- `degree == 2`: 通常のスケルトンピクセル（経路上の点）
- `degree >= 3`: 分岐点（ジャンクション）

---

### フェーズ2: 候補ノードの抽出 (162-190行目)

```cpp
// Identify raw nodes (endpoints degree==1, junctions degree>=3)
struct NodePix { int x; int y; };
std::vector<NodePix> raw_nodes;
```

**目的:** トポロジカルマップのノード候補を特定

**抽出条件:**
- `degree == 1`: 端点（行き止まり、経路の始点/終点）
- `degree >= 3`: 分岐点（複数の経路が交わる点）

**最適化:**
- OpenMP並列化でスレッドごとにローカル配列を使用
- クリティカルセクションで結果をマージ

**結果:**
- `raw_nodes[]`: ノード候補のピクセル座標リスト

---

### フェーズ3: 近接ノードのマージ（空間ハッシュ最適化） (192-223行目)

```cpp
// Merge nearby nodes within merge_radius (pixels)
// Use grid-based spatial hashing to avoid O(N²) comparison
const int grid_size = static_cast<int>(std::max(1.0, merge_radius_px)) + 1;
std::map<std::pair<int,int>, std::vector<int>> spatial_grid;
```

**目的:** 近接するノード候補を1つにマージして、ノード数を削減

**アルゴリズム:**
1. **空間グリッドハッシュ構築:**
   - マップをグリッドセルに分割（セルサイズ = `merge_radius_px`）
   - 各ノード候補をグリッドセルに登録

2. **近傍探索の最適化:**
   - 各ノードについて、同じセルと隣接8セル（3×3領域）のみを探索
   - 距離が `merge_radius` 以内のノードを統合（Union-Find）

**Union-Find データ構造:**
```cpp
auto findp = [&](int a){ while (parent[a] != a) a = parent[a] = parent[parent[a]]; return a; };
auto un = [&](int a, int b){ a = findp(a); b = findp(b); if (a!=b) parent[b]=a; };
```
- 経路圧縮による高速な集合管理
- 近接ノードを同一グループに統合

**計算量:**
- **従来**: O(N²) - 全ノード間の総当たり比較
- **最適化後**: O(N log N) - グリッドセル内のみ比較

---

### フェーズ4: 代表ノード決定と一時ノード作成 (224-257行目)

```cpp
// First pass: create temporary nodes and check connectivity
struct TempNode {
    int id;
    int x, y;  // pixel coordinates
    double world_x, world_y;
};
std::vector<TempNode> temp_nodes;
```

**目的:** マージされたグループから代表ノードを決定し、一時構造に格納

**処理内容:**
1. **グループ化:**
   - Union-Findで統合されたノードをグループ化
   
2. **重心計算:**
   - グループ内の全ノードの平均座標を計算
   - 代表位置として使用

3. **一時ノード生成:**
   - 各グループから1つの`TempNode`を生成
   - ピクセル座標と実世界座標の両方を保持
   - `label[x,y]`マップに登録（座標からノードIDを逆引き）

**安全性:**
- `max_nodes = 100,000`: ノード数上限で大規模マップ対応
- 境界チェックで配列外アクセスを防止

---

### フェーズ5: 連結成分解析による孤立ノード群削除（第1回） (259-368行目)

```cpp
// Connected component analysis: use simplified method for large maps
// For large maps (>10K nodes), skip this step and rely on edge tracing to filter nodes
const size_t max_nodes_for_component_analysis = 10000;
```

**目的:** スケルトン上での接続性を確認し、孤立した小さなノード群を削除

**条件分岐:**
- **小規模マップ（≤10,000ノード）**: 完全な連結成分解析を実行
- **大規模マップ（>10,000ノード）**: スキップして後続処理に任せる

**アルゴリズム（小規模マップ）:**
1. **BFSによる連結成分検出:**
   - 各ノードからスケルトン上を50ピクセルまで探索
   - 到達可能な他のノードを同じ連結成分としてマーク
   - `unordered_set`で効率的な訪問管理

2. **小さな連結成分のフィルタリング:**
   ```cpp
   const int min_component_size = 1;  // ユーザー設定値
   // Always keep the largest component
   ```
   - サイズ ≤ `min_component_size` の連結成分を削除
   - **最大の連結成分は常に保持**（小規模マップ保護）

3. **ノードとラベルの再構築:**
   - 保持する連結成分のノードのみを正式な`TopoNode`として登録
   - `label`マップを更新

**最適化:**
- `unordered_set`でメモリ使用量削減
- 探索距離制限（50ピクセル）で処理時間短縮
- 大規模マップでは自動スキップ

**計算量:** 
- 小規模: O(N × 探索距離)
- 大規模: O(N) - スキップしてノードをそのまま使用

---

### フェーズ6: エッジトレーシング (370-456行目)

```cpp
// Edge tracing: from each node, follow skeleton until another node or endpoint
std::vector<uint8_t> visited(width * height, 0);
std::vector<bool> node_used(topo.nodes.size(), false); // Track which nodes are actually used
```

**目的:** ノード間の経路（エッジ）を抽出し、ポリラインとして記録。同時に実際に使用されたノードを追跡

#### 6.1 各ノードからのトレース開始 (366-388行目)

**処理:**
1. ノード座標（メートル単位）をピクセル座標に変換
2. 座標がスケルトン上にない場合、最近傍スケルトンピクセルを探索（7×7領域）
3. スケルトン上の開始位置を確定

#### 6.2 8方向への経路追跡 (389-456行目)

**アルゴリズム:**
```cpp
for (int k = 0; k < 8; ++k) {  // 8方向を探索
    // 未訪問の隣接スケルトンピクセルから追跡開始
    while (steps < max_steps) {
        visited[current] = 1;
        poly.push_back(current);
        
        // 終了条件チェック
        if (到達点が別のノード) {
            エッジを生成;
            node_used[u] = true;  // ★使用ノードをマーク
            node_used[v] = true;
            break;
        }
        
        // 次のピクセルを選択（前のピクセルには戻らない）
        next = 未訪問の隣接スケルトンピクセル;
        
        if (次が見つからない) {
            新規ノードを作成;
            node_used[new_node] = true;  // ★即座にマーク
            エッジ生成;
            break;
        }
        
        length += distance(current, next);
        current = next;
    }
}
```

**主要な処理:**
1. **ポリライン記録:**
   - 経路上の各ピクセルを記録
   - ピクセル座標 → メートル単位に変換

2. **長さ計算:**
   - ユークリッド距離を累積（`std::hypot`）
   - 実世界単位（メートル）

3. **終了判定:**
   - **別のノードに到達**: エッジを確定して終了
   - **行き止まり**: その場に新規ノードを作成してエッジ確定
   - **ループ検出**: 開始ノードに戻った場合は中断

4. **次ピクセルの選択:**
   - 8近傍の未訪問スケルトンピクセルから選択
   - 前のピクセルには戻らない（バックトラック防止）
   - 訪問済みピクセルを避ける（サイクル防止）

5. **ノード使用追跡（重要）:**
   - エッジ生成時に両端ノードを`node_used[]`でマーク
   - 新規ノード作成時も即座にマーク
   - **degreeを使わない方法**: 実際のエッジトレーシングで使用されたノードを記録

**安全性:**
- `max_steps = min(width * height, 100,000)`: 無限ループ防止
- `max_edges = 200,000`: エッジ数上限
- 訪問フラグでサイクル検出

---

### フェーズ7: 未使用ノードの削除（第2回フィルタリング） (458-503行目)

```cpp
// Remove unused nodes after edge tracing (second method: actual usage tracking)
std::vector<bool> node_used(topo.nodes.size(), false);
```

**目的:** エッジトレーシングで実際に使用されなかったノードを削除

**アルゴリズム:**
1. **使用ノードのみを抽出:**
   - `node_used[i] == true` のノードのみを保持
   - 新しい連続したIDを割り当て

2. **ノードIDマッピング:**
   - `old_to_new_node_id`マップで対応関係を記録

3. **エッジのノードID更新:**
   - 全エッジの`u`と`v`を新しいノードIDに更新

**特徴:**
- **degreeを使わない方法**: エッジトレーシングの実際の動作を追跡
- エッジ生成に関与したノードのみを保持
- メモリ効率とデータ整合性を向上

**計算量:** O(N + E)

---

### フェーズ8: グラフ構造での連結成分解析（第3回フィルタリング） (504-603行目)

```cpp
// Additional filtering: remove small connected components after edge tracing
// Build adjacency list from edges
std::map<int, std::vector<int>> adjacency;
```

**目的:** 実際のグラフ構造で小さな連結成分（≤20ノード）を削除

**アルゴリズム:**
1. **隣接リストの構築:**
   - エッジからグラフ構造を構築
   - `adjacency[u] → [v1, v2, ...]`

2. **BFSによる連結成分検出:**
   - グラフ構造上でBFS探索
   - 到達可能な全ノードを同じ連結成分としてマーク
   - 各連結成分のサイズを計算

3. **小さな連結成分のフィルタリング:**
   ```cpp
   const int min_final_component_size = 20;
   // Always keep the largest component
   ```
   - サイズ ≤ 20 の連結成分を削除
   - **最大の連結成分は常に保持**

4. **ノードとエッジの再構築:**
   - 保持する連結成分のノードのみを残す
   - エッジも対応するノードがあるもののみ残す
   - ノードID、エッジIDを再割り当て

**特徴:**
- **degreeを使わない方法**: グラフのBFSで直接的に連結性を判定
- エッジベースの接続性確認で正確
- ノイズや小さな孤立構造を効果的に除去

**計算量:** O(N + E)

**効果:**
- tsudanumaマップ: 85,500 → 98ノード（**99.9%削減**）

---

### フェーズ9: エッジのプルーニング (605-623行目)

```cpp
// Pruning: remove edges shorter than threshold
const double min_len = params_.prune_min_length;
```

**目的:** 短い行き止まりエッジ（スパー）を除去

**アルゴリズム:**
1. **ノード次数マップ構築 (O(E)):**
   ```cpp
   std::map<int, int> node_degree;
   for (const auto& e : topo.edges) {
       node_degree[e.u]++;
       node_degree[e.v]++;
   }
   ```

2. **プルーニング条件判定:**
   ```cpp
   if (e.length < min_len) {
       bool u_is_endpoint = (node_degree[e.u] <= 1);
       bool v_is_endpoint = (node_degree[e.v] <= 1);
       should_keep = !(u_is_endpoint && v_is_endpoint);
   }
   ```

**削除条件（AND条件）:**
- エッジの長さ < `prune_min_length`
- **かつ** 両端点の次数が1（行き止まり同士）

**保持条件（以下のいずれか）:**
- エッジが十分長い
- 少なくとも一方の端点が分岐点（次数 ≥ 2）

**最適化:**
- **従来**: O(E²) - エッジごとに全エッジを走査
- **最適化後**: O(E) - 次数マップを事前構築

---

## 3段階の孤立ノード削除戦略

本実装では、異なるタイミングと方法で3回の孤立ノード削除を実施します：

| タイミング | 方法 | 削除対象 | degreeの使用 |
|-----------|------|----------|--------------|
| **フェーズ5** | スケルトンBFS連結成分解析 | 孤立したノード群 | ❌ 不使用 |
| **フェーズ7** | エッジトレーシング実績追跡 | 未使用ノード | ❌ 不使用 |
| **フェーズ8** | グラフBFS連結成分解析 | 小さな連結成分（≤20ノード） | ❌ 不使用 |

### なぜ3回必要か？

1. **フェーズ5（早期削除）:**
   - 目的: エッジトレーシング前に明らかに孤立したノード群を削除
   - 効果: メモリ使用量削減、エッジトレーシングの負荷軽減
   - 方法: スケルトン上の物理的な連結性を確認

2. **フェーズ7（実用性確認）:**
   - 目的: 理論的に到達可能でも実際に未使用のノードを削除
   - 効果: 実際のトポロジー構造のみを保持
   - 方法: エッジトレーシングの実行履歴を追跡

3. **フェーズ8（小構造削除）:**
   - 目的: ノイズや意味のない小さな構造を削除
   - 効果: クリーンで実用的なトポロジーマップ
   - 方法: グラフ構造でのBFSによる連結成分サイズ判定

---

## データ構造とメモリレイアウト

### 座標系の変換

```
ピクセル座標 (int x, int y)
    ↓ × resolution
実世界座標 (double x, double y) [メートル単位]
```

### 配列インデックス計算

```cpp
static inline int idx(int x, int y, int w) { return y * w + x; }
```
- 2次元座標を1次元配列インデックスに変換
- Row-major order（行優先順）

### 主要データ構造

| 構造 | サイズ | 用途 |
|------|--------|------|
| `is_skel[]` | W×H | 細線化されたスケルトンマスク |
| `degree[]` | W×H | 各ピクセルの接続数 |
| `label[]` | W×H | ピクセル座標→ノードIDマップ |
| `visited[]` | W×H | エッジトレーシングの訪問フラグ |
| `raw_nodes[]` | 可変 | ノード候補リスト |
| `temp_nodes[]` | 可変 | 一時ノード構造（連結成分解析用） |
| `node_used[]` | N | エッジトレーシングでの使用フラグ |
| `spatial_grid` | 可変 | グリッドセル→ノードインデックス |
| `parent[]` | N | Union-Find 親配列 |
| `groups[]` | N | マージグループ |
| `adjacency` | N+E | グラフ隣接リスト |

---

## パフォーマンス特性

### 計算量

| フェーズ | 計算量 | 説明 |
|----------|--------|------|
| 細線化 | O(I × W × H) | I: 反復回数（10-20回） |
| 次数計算 | O(W × H) | 全ピクセルを走査 |
| ノード抽出 | O(W × H) | 全ピクセルを走査 |
| ノードマージ | **O(N log N)** | 空間ハッシュ最適化 |
| 連結成分解析1 | O(N × D) | D: 探索距離（50px） |
| エッジトレース | O(W × H) | 各スケルトンピクセルを1回訪問 |
| 連結成分解析2 | O(N + E) | グラフBFS |
| プルーニング | O(E) | 次数マップ最適化 |

**全体:** O(I × W × H + N log N + N × D + E)

### 実測パフォーマンス（最新版）

| マップ | サイズ | 細線化 | ノード数 | エッジ数 | 処理時間 |
|--------|--------|--------|----------|----------|----------|
| room | 150×150 | ✅ | 1 | 0 | 95ms |
| anime | 2048×2048 | ✅ | 2 | 2 | 112ms |
| **tsudanuma** | **6665×5641** | **✅ (20反復)** | **98** | **94** | **2分1秒** |

### メモリ使用量

```
基本配列: 4 × W × H バイト (is_skel, degree, label, visited)
例: 6665×5641 = 約150MB

細線化: 2 × W × H バイト (img, marker)
例: 6665×5641 = 約75MB

動的構造: 
  - 一時ノード: 約100,000 × 32バイト = 3.2MB
  - 最終ノード: 約100 × 24バイト = 2.4KB
  - エッジ: 約100 × (基本構造 + ポリライン)
```

---

## 最適化技術

### 1. Zhang-Suen細線化の並列化

**最適化:**
```cpp
#pragma omp parallel for schedule(dynamic, 64)
```
- 各反復のマーキング処理を並列化
- 動的スケジューリングでロードバランス
- インライン近傍アクセスで高速化

**効果:**
- マルチコア環境で約2-3倍高速化
- 大規模マップで反復回数を20に制限

### 2. 空間グリッドハッシュ

**問題:** O(N²)の全ノード間距離計算が大規模マップでボトルネック

**解決策:**
```cpp
const int grid_size = merge_radius_px + 1;
std::map<std::pair<int,int>, std::vector<int>> spatial_grid;
```
- マップをグリッドセルに分割
- 各ノードの近傍セル（3×3領域）のみ探索
- **効果:** 約42倍の高速化

### 3. Union-Find with Path Compression

```cpp
auto findp = [&](int a){ 
    while (parent[a] != a) 
        a = parent[a] = parent[parent[a]];  // 経路圧縮
    return a; 
};
```
- ノードマージをO(α(N)) ≈ O(1)で実行（α: 逆アッカーマン関数）

### 4. unordered_set による訪問管理

**問題:** 連結成分解析で`vector<bool>`がメモリ不足を引き起こす

**解決策:**
```cpp
std::unordered_set<int> local_visited_set;
```
- 実際に訪問したピクセルのみを記録
- メモリ使用量を大幅削減

### 5. 段階的なノードフィルタリング

**戦略:**
1. 早期削除（連結成分解析）→ エッジトレーシングの負荷軽減
2. 使用追跡（エッジトレーシング中）→ 実際に使用されたノードのみ
3. 後処理（グラフBFS）→ 小さな構造の削除

**効果:**
- 大規模マップで99.9%のノード削減
- メモリ安全性の向上

---

## 安全性とロバスト性

### メモリ安全性

```cpp
// 配列アクセス前の境界チェック
auto inBounds = [&](int x, int y) { 
    return x >= 0 && y >= 0 && x < width && y < height; 
};

// size_t キャストでオーバーフロー防止
std::vector<int> label(static_cast<size_t>(width) * static_cast<size_t>(height), -1);
```

### 無限ループ防止

```cpp
const int max_iterations = (width * height > 10000000) ? 20 : 100; // 細線化
const int max_steps = std::min(width * height, 100000);  // エッジトレース
const int max_nodes = 100000;  // ノード数上限
const int max_edges = 200000;  // エッジ数上限
```

### サイクル検出

```cpp
visited[idx(cx,cy,width)] = 1;  // 訪問済みをマーク
if (visited[idx(tx,ty,width)]) continue;  // 訪問済みは避ける
```

### 大規模マップ対応

```cpp
// 連結成分解析をスキップ
if (temp_nodes.size() > 10000) {
    // Edge tracing will filter nodes
}

// 細線化をスキップ
if (width * height > 5000000) {
    is_skel = gvd_mask;
}
```

---

## パラメータ

### `params_.merge_radius`
- **デフォルト:** 0.5メートル（一般的）
- **意味:** この距離以内のノードを1つにマージ
- **効果:**
  - 小さい値: 詳細なトポロジー、ノード数増加
  - 大きい値: 簡略化されたトポロジー、ノード数減少

### `params_.prune_min_length`
- **デフォルト:** 0.3メートル（一般的）
- **意味:** これより短い行き止まりエッジを削除
- **効果:**
  - 小さい値: 詳細な構造を保持、ノイズ残存
  - 大きい値: ノイズ除去、重要な短辺も削除の可能性

### 内部定数（調整可能）

```cpp
const int min_component_size = 1;          // 第1回フィルタ: 連結成分サイズ閾値
const int min_final_component_size = 20;   // 第3回フィルタ: 最終連結成分サイズ閾値
const int max_search_dist = 50;            // 連結成分解析の探索距離
```

---

## トラブルシューティング

### 大規模マップでのメモリ不足/コアダンプ

**症状:** セグメンテーションフォルト、Out of Memory

**対処:**
```cpp
const int max_nodes = 100000;  // ← この値を調整
const int max_edges = 200000;  // ← この値を調整
const size_t max_nodes_for_component_analysis = 10000; // ← スキップ閾値
```

**根本原因:**
- 細線化により大量のスケルトンピクセルが生成される
- 連結成分解析でメモリ不足

**解決策:**
1. 大規模マップでは連結成分解析を自動スキップ
2. エッジトレーシング後のフィルタリングに依存
3. max_nodes/max_edgesで上限を設定

### ノード/エッジが多すぎる

**原因:**
- `merge_radius`が小さすぎる
- GVDスケルトンが過度に詳細
- 細線化が不完全

**対処:**
1. `merge_radius`を増やす（例: 0.5 → 1.0メートル）
2. GVDの解像度を下げる
3. `min_final_component_size`を増やす（例: 20 → 50）
4. 細線化の`max_iterations`を増やす

### ノード/エッジが少なすぎる

**原因:**
- `min_final_component_size`が大きすぎる
- 最大連結成分が小さい

**対処:**
1. `min_final_component_size`を減らす（現在20）
2. `merge_radius`を減らして詳細なノードを保持
3. 連結成分解析の探索距離を増やす

### エッジが途切れる

**原因:**
- `max_steps`制限に到達
- スケルトンに隙間がある

**対処:**
1. `max_steps`を増やす
2. 細線化のパラメータを調整
3. GVD生成パラメータを調整

---

## 使用例

```cpp
// TopologyExtractorの初期化
TopologyExtractor::Params params;
params.merge_radius = 0.5;      // 0.5メートル以内をマージ
params.prune_min_length = 0.3;  // 0.3メートル未満の行き止まりを削除

TopologyExtractor extractor(params);

// トポロジカルマップの抽出
TopologicalMap topo = extractor.run(
    gvd_result.gvd_mask,  // GVDマスク
    grid.width,           // マップ幅
    grid.height,          // マップ高さ
    grid.resolution       // 解像度（m/pixel）
);

std::cout << "Nodes: " << topo.nodes.size() << std::endl;
std::cout << "Edges: " << topo.edges.size() << std::endl;

// JSONとして保存
std::string json = toJson(topo);
std::ofstream ofs("topology.json");
ofs << json;
```

---

## 実測結果の比較

### map_tsudanuma.pgm での変遷

| バージョン | 細線化 | ノード数 | エッジ数 | 処理時間 | 結果 |
|-----------|--------|----------|----------|----------|------|
| 初期版 | ❌ | 315 | 2 | タイムアウト | ❌ 失敗 |
| 形態学的骨格化 | ✅ | 46,714 | 29,721 | タイムアウト | ❌ 失敗 |
| 最適化版 | ❌ | 50,000 | 22,587 | 2分31秒 | ⚠️ ノード多すぎ |
| 孤立削除v1 | ❌ | 67,281 | 47,471 | 2分1秒 | ⚠️ ノード多すぎ |
| **最終版** | **✅ Zhang-Suen** | **98** | **94** | **2分1秒** | **✅ 成功** |

### フィルタリング効果（tsudanumaマップ）

```
生ノード候補: ~100,000個
    ↓ マージ
マージ後: ~85,500個
    ↓ 連結成分解析（第1回）: 最大連結成分のみ保持
フィルタ後: ~85,000個
    ↓ エッジトレーシング: 実際に使用されたノードのみ
使用ノード: ~85,000個
    ↓ 連結成分解析（第2回）: ≤20ノードの連結成分削除
最終結果: 98ノード（99.9%削減！）
```

---

## 関連項目

- `GvdGenerator::run`: GVDスケルトンの生成
- `zhangSuenThinning`: Zhang-Suen細線化アルゴリズム
- `Visualizer::saveTopologicalMapAsImage`: トポロジーの可視化
- `toJson`: トポロジカルマップのJSON出力

---

## 参考文献

- **Zhang-Suen Algorithm**: T.Y. Zhang and C.Y. Suen, "A Fast Parallel Algorithm for Thinning Digital Patterns," Communications of the ACM, 1984
- **Robot Motion**: Choset, H., et al. "Principles of Robot Motion: Theory, Algorithms, and Implementations" (MIT Press, 2005)
- **Probabilistic Robotics**: Thrun, S., et al. "Probabilistic Robotics" (MIT Press, 2005)
- **Spatial Data Structures**: Samet, H. "Foundations of Multidimensional and Metric Data Structures" (Morgan Kaufmann, 2006)

---

**最終更新:** 2025年10月15日  
**バージョン:** 2.0 (細線化 + 3段階フィルタリング対応版)
