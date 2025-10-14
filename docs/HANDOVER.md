# General Voronoi Diagram プロジェクト 引き継ぎ説明書

## プロジェクト概要

このプロジェクトは、占有格子地図から一般化ボロノイ図（GVD）を使用してトポロジカルマップを生成するC++ライブラリです。ROS互換の設計となっており、ロボットナビゲーションでの使用を想定しています。

## プロジェクト構造

```
General_Voronoi_Diagram/
├── CMakeLists.txt              # CMakeビルド設定
├── README.md                   # 英語版README
├── README_JP.md               # 日本語版README
├── HANDOVER.md                # この引き継ぎ説明書
├── plan.txt                   # 開発計画書
├── LICENSE                    # ライセンスファイル
│
├── include/gvd_topo/          # ヘッダーファイル
│   ├── gvd_topo.hpp          # メインヘッダー（全機能のエントリーポイント）
│   ├── core/                 # コア機能
│   │   ├── OccupancyGrid.hpp
│   │   ├── GvdGenerator.hpp
│   │   ├── TopologyExtractor.hpp
│   │   └── Visualizer.hpp
│   ├── io/                   # 入出力機能
│   │   └── YamlLoader.hpp
│   ├── utils/                # ユーティリティ
│   │   ├── Timer.hpp
│   │   └── ConfigManager.hpp
│   ├── cli/                  # CLI機能
│   │   └── CliApplication.hpp
│   ├── ros_adapters.hpp      # ROS互換アダプタ
│   └── parameters.hpp        # パラメータ定義
│
├── src/                      # ソースファイル
│   ├── main.cpp             # ライブラリエントリーポイント
│   ├── core/                # コア機能実装
│   │   ├── OccupancyGrid.cpp
│   │   ├── GvdGenerator.cpp
│   │   ├── TopologyExtractor.cpp
│   │   └── Visualizer.cpp
│   ├── io/                  # 入出力機能実装
│   │   └── YamlLoader.cpp
│   ├── utils/               # ユーティリティ実装
│   │   └── ConfigManager.cpp
│   ├── cli/                 # CLI機能実装
│   │   ├── cli_main.cpp
│   │   └── CliApplication.cpp
│   ├── ros_adapters.cpp     # ROS互換アダプタ実装
│   └── parameters.cpp       # パラメータ実装
│
├── tools/                   # ツール類
│   ├── create_test_maps.cpp
│   ├── regression_test.cpp
│   ├── optimized_new_map.cpp
│   └── debug_tools/         # デバッグツール
│       ├── debug_topology.cpp
│       ├── debug_visualization.cpp
│       └── debug_new_map.cpp
│
├── tests/                   # テストファイル
│   ├── unit/               # 単体テスト
│   │   └── test_ros_adapters.cpp
│   ├── integration/        # 統合テスト（将来実装）
│   └── test_data/          # テストデータ（将来実装）
│
├── data/                   # データファイル
│   ├── example_config.yaml # サンプル設定ファイル
│   ├── *.pgm              # サンプルマップ画像
│   ├── *.yaml             # マップ設定ファイル
│   └── *.json             # 生成されたトポロジカルマップ
│
└── build/                  # ビルド出力（生成される）
    ├── bin/               # 実行ファイル
    └── lib/               # ライブラリファイル
```

## 各ファイルの役割

### コア機能（src/core/）

#### OccupancyGrid.cpp
- **役割**: 占有格子地図の管理と操作
- **主要機能**:
  - 画像ファイルからの読み込み（PNG/PGM）
  - ランダムマップの生成（ベンチマーク用）
  - 占有判定とデータ変換
  - デバッグ用PGMファイル出力

#### GvdGenerator.cpp
- **役割**: 一般化ボロノイ図（GVD）の生成
- **主要機能**:
  - ユークリッド距離変換（EDT）の実行
  - リッジ検出によるGVD抽出
  - OpenCVを使用した高速処理
  - バイナリマスクとポリラインの出力

#### TopologyExtractor.cpp
- **役割**: GVDからトポロジカルマップの抽出
- **主要機能**:
  - ノード（分岐点・終端点）の検出
  - エッジ（経路）の抽出
  - 不要な枝の刈り取り
  - 近接ノードのマージ

#### Visualizer.cpp
- **役割**: トポロジカルマップの可視化
- **主要機能**:
  - トポロジカルマップの画像生成
  - ノードとエッジの描画
  - 設定可能な可視化オプション

### 入出力機能（src/io/）

#### YamlLoader.cpp
- **役割**: YAMLファイルの読み込みと設定管理
- **主要機能**:
  - マップ情報のYAML読み込み
  - 設定ファイルの読み込み・保存
  - 相対パスの絶対パス変換
  - 設定の検証

### ユーティリティ（src/utils/）

#### ConfigManager.cpp
- **役割**: 設定の統合管理
- **主要機能**:
  - コマンドライン引数の解析
  - 設定ファイルとの統合
  - 設定の検証
  - ヘルプメッセージの表示

### CLI機能（src/cli/）

#### cli_main.cpp
- **役割**: コマンドラインアプリケーションのエントリーポイント
- **主要機能**:
  - 設定管理の初期化
  - アプリケーションの実行
  - エラーハンドリング

#### CliApplication.cpp
- **役割**: CLIアプリケーションのメインロジック
- **主要機能**:
  - 処理パイプラインの実行
  - タイミング測定
  - 出力ファイルの生成
  - 統計情報の表示

### ツール類（tools/）

#### create_test_maps.cpp
- **役割**: テスト用マップの生成
- **生成するマップ**:
  - 廊下マップ
  - 部屋マップ
  - T字路マップ

#### regression_test.cpp
- **役割**: 回帰テストの実行
- **機能**:
  - 既知のマップでのテスト実行
  - 結果の比較
  - ゴールデンテストの生成

#### debug_*.cpp
- **役割**: デバッグ用ツール
- **機能**:
  - 特定の処理段階のデバッグ
  - 中間結果の可視化
  - パフォーマンス測定

## ビルド方法

### 前提条件
- C++17対応コンパイラ
- CMake 3.16以上
- OpenCV（オプション、推奨）
- OpenMP（オプション）

### ビルド手順

```bash
# ビルドディレクトリの作成と設定
cmake -B build -S .

# ビルド実行
cmake --build build -j

# または
make -C build
```

### ビルドオプション

```bash
# OpenCV無効でビルド
cmake -B build -S . -DWITH_OPENCV=OFF

# デバッグツール無効でビルド
cmake -B build -S . -DBUILD_TOOLS=OFF

# ネイティブ最適化無効でビルド
cmake -B build -S . -DUSE_NATIVE_OPTIMIZATIONS=OFF
```

## 実行方法

### 基本的な使用方法

```bash
# ヘルプの表示
./build/bin/gvd_topo_cli --help

# 画像ファイルからトポロジカルマップ生成
./build/bin/gvd_topo_cli --input data/room.pgm --resolution 0.05 --out-map room_topo.json

# YAML設定ファイルを使用
./build/bin/gvd_topo_cli --config data/example_config.yaml

# ベンチマークモード
./build/bin/gvd_topo_cli --bench-w 100 --bench-h 100 --out-map benchmark.json
```

### 設定ファイルの使用

`data/example_config.yaml`を参考に設定ファイルを作成できます：

```yaml
input:
  image_file: "room.pgm"
  resolution: 0.05
  occupancy_threshold: 50

output:
  map_file: "output_topo.json"
  gvd_image: "output_gvd.png"
  topo_image: "output_topo.png"

benchmark:
  enabled: false
  width: 100
  height: 100
  occupancy_ratio: 0.3
  seed: 12345

processing:
  distance_epsilon: 1e-6
  use_opencv: true
  prune_min_length: 0.5
  merge_radius: 0.2
  max_trace_steps: 100000
```

### ライブラリとしての使用

```cpp
#include "gvd_topo/gvd_topo.hpp"

// 占有格子地図の読み込み
auto grid = gvd_topo::OccupancyGrid::loadFromImage("map.pgm", 0.05, 50);

// GVD生成
gvd_topo::GvdGenerator gvd;
auto gvd_result = gvd.run(grid);

// トポロジー抽出
gvd_topo::TopologyExtractor topo;
auto topo_map = topo.run(gvd_result.gvd_mask, gvd_result.width, gvd_result.height, 0.05);

// 結果の保存
std::ofstream ofs("output.json");
ofs << gvd_topo::toJson(topo_map) << std::endl;
```

## 主要なクラスとAPI

### OccupancyGrid
```cpp
class OccupancyGrid {
public:
    // コンストラクタ
    OccupancyGrid(int width, int height, double resolution);
    
    // 静的メソッド
    static OccupancyGrid loadFromImage(const std::string& path, double resolution, int threshold);
    static OccupancyGrid randomMap(int width, int height, double resolution, double occupancy_ratio, unsigned seed);
    
    // メンバ変数
    int width, height;
    double resolution;
    std::vector<int8_t> data;
    // ...
};
```

### GvdGenerator
```cpp
class GvdGenerator {
public:
    GvdResult run(const OccupancyGrid& grid);
};

struct GvdResult {
    std::vector<bool> gvd_mask;
    int width, height;
    // ...
};
```

### TopologyExtractor
```cpp
class TopologyExtractor {
public:
    TopologicalMap run(const std::vector<bool>& gvd_mask, int width, int height, double resolution);
};

struct TopologicalMap {
    std::vector<TopoNode> nodes;
    std::vector<TopoEdge> edges;
    // ...
};
```

## デバッグとテスト

### デバッグツールの使用

```bash
# トポロジー抽出のデバッグ
./build/bin/debug_topology

# 可視化のデバッグ
./build/bin/debug_visualization

# 新しいマップのデバッグ
./build/bin/debug_new_map
```

### テストの実行

```bash
# 回帰テスト
./build/bin/regression_test

# ROSアダプタのテスト
./build/bin/test_ros_adapters

# テストマップの生成
./build/bin/create_test_maps
```

## トラブルシューティング

### よくある問題

1. **ビルドエラー**
   - OpenCVが見つからない場合：`-DWITH_OPENCV=OFF`でビルド
   - コンパイラエラーの場合：C++17対応コンパイラを使用

2. **実行時エラー**
   - ファイルが見つからない：パスを確認
   - メモリ不足：マップサイズを小さくする

3. **設定ファイルエラー**
   - YAML構文エラー：インデントを確認
   - パスエラー：相対パスを確認

### ログとデバッグ

- タイミング情報は標準出力に表示
- エラーメッセージは標準エラーに表示
- デバッグツールで中間結果を確認可能

## 今後の開発

### 計画されている機能

1. **テストの充実**
   - 単体テストの追加
   - 統合テストの実装
   - パフォーマンステスト

2. **機能拡張**
   - より多くの画像形式対応
   - 並列処理の最適化
   - ROS2対応

3. **ドキュメント**
   - APIドキュメントの生成
   - 使用例の追加
   - チュートリアルの作成

### 開発ガイドライン

1. **コードスタイル**
   - クラス名：PascalCase
   - 関数名：camelCase
   - 定数名：UPPER_CASE
   - ファイル名：snake_case

2. **コミットメッセージ**
   - 機能追加：`feat: 機能の説明`
   - バグ修正：`fix: 修正内容`
   - リファクタリング：`refactor: 変更内容`

3. **テスト**
   - 新機能追加時はテストも追加
   - 既存テストが通ることを確認
   - パフォーマンステストも考慮

## 連絡先・サポート

- プロジェクトの詳細は`README.md`を参照
- 開発計画は`plan.txt`を参照
- 問題が発生した場合は、ログとエラーメッセージを確認

---

**最終更新**: 2024年9月25日
**バージョン**: 1.0.0
**ライセンス**: プロジェクトのLICENSEファイルを参照

