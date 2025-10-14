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

### フェーズ1: スケルトンピクセルの次数計算 (23-42行目)

```cpp
// Compute degrees for each skeleton pixel (8-neighborhood)
std::vector<uint8_t> degree(width * height, 0);
const int dx8[8] = { -1, 0, 1, -1, 1, -1, 0, 1 };
const int dy8[8] = { -1,-1,-1,  0, 0,  1, 1, 1 };
```

**目的:** 各スケルトンピクセルの8近傍における接続数（次数）を計算

**処理内容:**
- 各スケルトンピクセルについて、8方向の隣接ピクセルを確認
- スケルトンピクセルの隣接数をカウント → `degree[]` に格納
- OpenMP並列化により高速処理（`#pragma omp parallel for`）

**次数の意味:**
- `degree == 1`: 端点（エンドポイント）
- `degree == 2`: 通常のスケルトンピクセル（経路上の点）
- `degree >= 3`: 分岐点（ジャンクション）

---

### フェーズ2: 候補ノードの抽出 (44-72行目)

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

### フェーズ3: 近接ノードのマージ（空間ハッシュ最適化） (74-111行目)

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

### フェーズ4: 代表ノードの決定とラベル付け (112-129行目)

```cpp
// Compute representatives and average positions
std::vector<std::vector<int>> groups(raw_nodes.size());
for (size_t i = 0; i < raw_nodes.size(); ++i) 
    groups[findp(i)].push_back(i);
```

**目的:** マージされたグループから代表ノードを決定

**処理内容:**
1. **グループ化:**
   - Union-Findで統合されたノードをグループ化
   
2. **重心計算:**
   - グループ内の全ノードの平均座標を計算
   - 代表位置として使用

3. **ノード生成:**
   - 各グループから1つの`TopoNode`を生成
   - ピクセル座標 → メートル単位に変換（`× resolution`）
   - `label[x,y]`マップに登録（座標からノードIDを逆引き）

**安全性:**
- `max_nodes = 50,000`: ノード数上限で大規模マップ対応
- 境界チェックで配列外アクセスを防止

---

### フェーズ5: エッジトレーシング (131-225行目)

```cpp
// Edge tracing: from each node, follow skeleton until another node or endpoint
std::vector<uint8_t> visited(width * height, 0);
```

**目的:** ノード間の経路（エッジ）を抽出し、ポリラインとして記録

#### 5.1 各ノードからのトレース開始 (141-164行目)

**処理:**
1. ノード座標（メートル単位）をピクセル座標に変換
2. 座標がスケルトン上にない場合、最近傍スケルトンピクセルを探索（7×7領域）
3. スケルトン上の開始位置を確定

#### 5.2 8方向への経路追跡 (165-224行目)

**アルゴリズム:**
```cpp
for (int k = 0; k < 8; ++k) {  // 8方向を探索
    // 未訪問の隣接スケルトンピクセルから追跡開始
    while (steps < max_steps) {
        visited[current] = 1;
        poly.push_back(current);
        
        // 終了条件チェック
        if (到達点が別のノード) {
            エッジを生成して終了;
        }
        
        // 次のピクセルを選択（前のピクセルには戻らない）
        next = 未訪問の隣接スケルトンピクセル;
        
        if (次が見つからない) {
            新規ノードを作成してエッジ生成;
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

**安全性:**
- `max_steps = min(width * height, 100,000)`: 無限ループ防止
- `max_edges = 100,000`: エッジ数上限
- 訪問フラグでサイクル検出

**生成されるエッジデータ:**
```cpp
TopoEdge {
    id: エッジID,
    u: 始点ノードID,
    v: 終点ノードID,
    length: 実世界距離（メートル）,
    polyline: [(x1,y1), (x2,y2), ...] // メートル単位
}
```

---

### フェーズ6: エッジのプルーニング (227-251行目)

```cpp
// Pruning: remove edges shorter than threshold
const double min_len = params_.prune_min_length;
```

**目的:** ノイズや短いスパー（行き止まり）を除去

**アルゴリズム:**
1. **ノード次数マップ構築 (O(E)):**
   ```cpp
   std::map<int, int> node_degree;
   for (const auto& e : topo.edges) {
       node_degree[e.u]++;
       node_degree[e.v]++;
   }
   ```
   - 各ノードが接続するエッジ数をカウント

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
| `is_skel[]` | W×H | スケルトンマスク |
| `degree[]` | W×H | 各ピクセルの接続数 |
| `label[]` | W×H | ピクセル座標→ノードIDマップ |
| `visited[]` | W×H | エッジトレーシングの訪問フラグ |
| `raw_nodes[]` | 可変 | ノード候補リスト |
| `spatial_grid` | 可変 | グリッドセル→ノードインデックス |
| `parent[]` | N | Union-Find 親配列 |
| `groups[]` | N | マージグループ |

---

## パフォーマンス特性

### 計算量

| フェーズ | 計算量 | 説明 |
|----------|--------|------|
| 次数計算 | O(W×H) | 全ピクセルを走査 |
| ノード抽出 | O(W×H) | 全ピクセルを走査 |
| ノードマージ | **O(N log N)** | 空間ハッシュ最適化 |
| エッジトレース | O(W×H) | 各スケルトンピクセルを1回訪問 |
| プルーニング | **O(E)** | 次数マップ最適化 |

**全体:** O(W×H + N log N + E)

### 実測パフォーマンス

| マップ | サイズ | ノード数 | エッジ数 | 処理時間 |
|--------|--------|----------|----------|----------|
| room | 150×150 | 8 | 3 | 3.9ms |
| anime | 2048×2048 | 21 | 16 | 38ms |
| tsudanuma | 6665×5641 | 50,000 | 22,640 | **2.7秒** |

### メモリ使用量

```
基本配列: 4 × W × H バイト (is_skel, degree, label, visited)
例: 6665×5641 = 約150MB

動的構造: 
  - ノード: 約50,000 × 24バイト = 1.2MB
  - エッジ: 約22,640 × (基本構造 + ポリライン)
```

---

## 最適化技術

### 1. 空間グリッドハッシュ

**問題:** O(N²)の全ノード間距離計算が大規模マップでボトルネック

**解決策:**
```cpp
const int grid_size = merge_radius_px + 1;
std::map<std::pair<int,int>, std::vector<int>> spatial_grid;
```
- マップをグリッドセルに分割
- 各ノードの近傍セル（3×3領域）のみ探索
- **効果:** 約42倍の高速化（tsudanumaマップ）

### 2. Union-Find with Path Compression

```cpp
auto findp = [&](int a){ 
    while (parent[a] != a) 
        a = parent[a] = parent[parent[a]];  // 経路圧縮
    return a; 
};
```
- ノードマージをO(α(N)) ≈ O(1)で実行（α: 逆アッカーマン関数）

### 3. 次数マップの事前構築

**問題:** プルーニング時のO(E²)計算

**解決策:**
```cpp
std::map<int, int> node_degree;
for (const auto& e : topo.edges) {
    node_degree[e.u]++;
    node_degree[e.v]++;
}
```
- 一度だけ全エッジを走査してノード次数を記録
- プルーニング判定がO(1)に

### 4. OpenMP並列化

```cpp
#pragma omp parallel for schedule(static)
```
- 次数計算とノード抽出をマルチスレッド化
- 静的スケジューリングでキャッシュ効率向上

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
const int max_steps = std::min(width * height, 100000);
const int max_nodes = 50000;
const int max_edges = 100000;
```

### サイクル検出

```cpp
visited[idx(cx,cy,width)] = 1;  // 訪問済みをマーク
if (visited[idx(tx,ty,width)]) continue;  // 訪問済みは避ける
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

---

## トラブルシューティング

### 大規模マップでのメモリ不足

**症状:** セグメンテーションフォルト、Out of Memory

**対処:**
```cpp
const int max_nodes = 50000;  // ← この値を調整
const int max_edges = 100000; // ← この値を調整
```

### ノード/エッジが多すぎる

**原因:**
- `merge_radius`が小さすぎる
- GVDスケルトンが過度に詳細

**対処:**
1. `merge_radius`を増やす（例: 0.5 → 1.0メートル）
2. GVDの解像度を下げる
3. より積極的なプルーニング（`prune_min_length`を増やす）

### エッジが途切れる

**原因:**
- `max_steps`制限に到達
- スケルトンに隙間がある

**対処:**
1. `max_steps`を増やす
2. GVD生成パラメータを調整

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

## 関連項目

- `GvdGenerator::run`: GVDスケルトンの生成
- `Visualizer::saveTopologicalMapAsImage`: トポロジーの可視化
- `toJson`: トポロジカルマップのJSON出力

---

## 参考文献

- Choset, H., et al. "Principles of Robot Motion: Theory, Algorithms, and Implementations" (MIT Press, 2005)
- Thrun, S., et al. "Probabilistic Robotics" (MIT Press, 2005)
- 空間分割法: Samet, H. "Foundations of Multidimensional and Metric Data Structures" (Morgan Kaufmann, 2006)

---

**最終更新:** 2025年10月14日  
**バージョン:** 1.0 (大規模マップ対応版)

