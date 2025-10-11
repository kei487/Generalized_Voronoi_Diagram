# 一般化ボロノイ図からのトポロジカルマップ生成

占有格子地図から一般化ボロノイ図（GVD）を使用してトポロジカルマップを生成するC++ライブラリです。

## 機能

- **占有格子地図処理**: PNG/PGM画像とROS互換の占有格子地図の読み込み
- **GVD生成**: ユークリッド距離変換とリッジ検出によるスケルトン抽出
- **トポロジ抽出**: GVDからノード（分岐点/終端点）とエッジの識別
- **ROS対応**: ROS統合用のメッセージアダプタとパラメータ構造体
- **性能最適化**: OpenMP並列化とネイティブ最適化
- **複数出力形式**: JSONトポロジカルマップとPNG可視化オーバレイ

## ビルド要件

- C++17コンパイラ
- CMake 3.16+
- OpenCV（オプション、画像I/Oと距離変換用）
- OpenMP（オプション、並列化用）

## クイックスタート

```bash
# ビルド
cmake -B build -S .
cmake --build build -j

# サンプルマップで実行
./build/bin/gvd_topo_cli --input data/room.pgm --resolution 0.05 \
  --out-map room_topo.json --out-gvd room_gvd.png
```

## 使用方法

### build cmd

```
  mkdir build && cd build && cmake .. && make -j$(nproc)
```
### コマンドラインインターフェース

```bash
./build/bin/gvd_topo_cli [オプション]

オプション:
  --input, -i <パス>        入力画像ファイル（PNG/PGM）
  --resolution <浮動小数>   マップ解像度（メートル/ピクセル、デフォルト: 0.05）
  --occ-thresh <整数>       占有判定閾値 0-100（デフォルト: 50）
  --out-map <パス>          出力トポロジカルマップJSON
  --out-gvd <パス>          出力GVDオーバレイPNG
  --bench-w <整数>          ベンチマークモード: 合成マップ幅
  --bench-h <整数>          ベンチマークモード: 合成マップ高さ
  --bench-occ <浮動小数>    ベンチマークモード: 障害物比率 0.0-1.0
```

### ライブラリ使用

```cpp
#include "gvd_topo/OccupancyGrid.hpp"
#include "gvd_topo/GvdGenerator.hpp"
#include "gvd_topo/TopologyExtractor.hpp"

// 占有格子地図の読み込み
auto grid = OccupancyGrid::loadFromImage("map.pgm", 0.05, 50);

// GVD生成
GvdGenerator gvd;
auto gvd_result = gvd.run(grid);

// トポロジ抽出
TopologyExtractor topo;
auto topological_map = topo.run(gvd_result.gvd_mask, 
                               gvd_result.width, gvd_result.height, 0.05);
```

## テストデータ

`data/`に合成テストマップが利用可能です：
- `room.pgm`: ドア付き部屋（150x150）
- `corridor.pgm`: 障害物付き廊下（200x100）
- `t_junction.pgm`: T字交差点レイアウト（200x200）

## 出力形式

トポロジカルマップはJSON形式で保存されます：

```json
{
  "nodes": [
    {"id": 0, "x": 3.75, "y": 3.75},
    {"id": 1, "x": 7.25, "y": 2.85}
  ],
  "edges": [
    {"id": 0, "u": 0, "v": 1, "length": 4.2, "polyline": [[3.75,3.75], [7.25,2.85]]}
  ]
}
```

## 性能

1000x1000合成マップでのベンチマーク結果：
- 読み込み+前処理: ~170ms
- EDT+GVD: ~110ms
- トポロジ: 複雑さにより変動

## ROS統合

ライブラリはROS互換のデータ構造を提供します：

```cpp
#include "gvd_topo/ros_adapters.hpp"
#include "gvd_topo/parameters.hpp"

// ROSメッセージ形式に変換
OccupancyGridMsg msg = toMsg(occupancy_grid);
OccupancyGrid grid = fromMsg(msg);

// 可視化マーカーの生成
MarkerArray markers = topologicalMapToMarkers(topological_map);

// パラメータ管理
ProcessingParameters params = getDefaultProcessingParameters();
```

## アルゴリズム概要

1. **占有格子地図前処理**: 二値化と形態学処理
2. **距離変換**: ユークリッド距離変換（EDT）による各セルから最近障害物までの距離計算
3. **GVD抽出**: 距離マップの局所最大値からスケルトン検出
4. **ノード抽出**: 分岐点（次数≥3）と終端点（次数=1）の識別
5. **エッジ抽出**: ノード間の経路トレースと長さ計算
6. **プルーニング**: 短い枝の除去

## パラメータ調整

主要パラメータ：
- `occupancy_threshold`: 占有判定閾値（デフォルト: 50）
- `prune_min_length`: 除去する最小枝長（デフォルト: 0.5m）
- `merge_radius`: 近傍ノード結合半径（デフォルト: 0.2m）

## ライセンス

詳細はLICENSEファイルを参照してください。
