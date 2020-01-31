# Pattern Mining
非著者による[safe pattern pruning](https://arxiv.org/abs/1602.04548)の実装

著者らの実装は[こちら](https://github.com/takeuchi-lab/SafePatternPruning)にあるが, 
グラフマイニングにおいてバグを見つけたため, [gBoost](https://github.com/rkwitt/gboost)を元に作り直した. 

itemset, sequence, graphデータセットに対するLASSOをShooting Algorithmで解く. 


**特徴:**
- 最適化において, 探索木の保存を行っている ([tree.hh](http://tree.phi-sci.com/)を使用)
- ただの頻出パターンマイニングでは探索木の保存を行っていない

このソースコードをベースに何かを実装すれば, 木の保存をしたい場合でも, したくない場合でも参考になるであろう. 

## ItemsetMining
アイテム集合からターゲットに影響を与える部分アイテム集合をマイニングする. 

### Dataの形式
```
+1 1 2 3 4 6
-1 2 3 4 6 7 8
-1 1 3 4 5 6
+1 3 5 6 7
```

- 各行がトランザクションを表す
- 1列目が回帰のターゲット (実数)
- 2列目以降の整数がアイテムを表す
- アイテムは昇順に並んでなければならない
- 各行のアイテムは重複してはいけない

## SequenceMining
系列データからターゲットに影響を与える部分系列をマイニングする. 

### Dataの形式
```
+1 1 2 2 1 2
-1 2 2 1 2 1 2
-1 1 2 1 1 1
+1 1 2 1 2
```

- 各行がトランザクションを表す
- 1列目が回帰のターゲット (実数)
- 2列目以降の整数が系列を表す

## GraphMining
グラフデータからターゲットに影響を与える部分グラフをマイニングする. 

[gBoost](https://github.com/rkwitt/gboost)をベースに書き換えている. 
### Dataの形式
```
t # 0 1
v 0 1
v 1 1
v 2 2
v 3 3
e 0 1 1
e 1 2 1
e 1 3 1

t # 0 -1
v 0 1
v 1 1
v 2 2
v 3 3
e 0 1 1
e 1 2 1
e 1 3 1
e 3 0 1
```

- tの行について
    - `t # 0 ターゲット(実数)`
    - #と0は冗長であるが, ベースとした手法の形式に従っているだけなため, 無視していただきたい
- vの行について
    - `v 頂点ID 頂点ラベル`
    - 各行が各頂点に対応している
- eの行について
    - `e 頂点ID 頂点ID 辺ラベル`
    - 各行が各辺に対応している

これはDFSコードという形式で, グラフを深さ優先探索することで作られる. 
詳細は[gSpan](https://dm.kaist.ac.kr/kse625/resources/Yan_2002.pdf)を参考にされたい. 

なお, DFSコードへの変換プログラムは[こちら](https://github.com/birdwatcherYT/GraphDatasetConverter)で公開中. 

## コンパイル&実行方法
### コンパイル
> make

### 実行
> ./run filename

- run: プログラム名
- filename: ファイル名
    - `./data/`からロードされる

### ハイパーパラメータについて
```
const int splitNum = 100;
const int maxloop = INT_MAX;
const double eps = 1e-6;
const int freq = 5;
const double R=pow(0.01, 1.0/(splitNum-1));
const int minsup=1;
const int maxpat=10;
```
`main.cpp`内に記されているため, 適当に書き換えていただきたい. 
- splitNum: 正則化係数λの更新回数. 線形モデルの係数が`w = 0`となるλ (λmaxと呼ぶ)から始まる. 
- maxloop: 最適化の最大ループ数
- eps: 収束条件. 相対双対ギャップがeps以下で終了する
- freq: Dynamic Screeningの実行間隔 (Safe Screeningの説明は[こちら](https://qiita.com/birdwatcher/items/6c3f86693f02762d05b9)にまとめた)
- R: λの減少係数. `λ = λ * R` によってλが変化していく. この例では, λmaxと0.01λmaxの間を対数スケールで100分割するようなRを使っている. 
- minsup: 最小サポート
- maxpat: 最大パターン長

### 関数について
実装されている重要な関数について軽く説明. 

- mining: 頻出パターンマイニングの実装 (木の保存なし)
    - 頻出パターンとそのサポート, 転置した計画行列(0/1行列)が得られる
    - __mining: 再帰実装の本体
- createRoot: 木の保存ありで, ルートノードを作る. 
- createChildren: 木の保存ありで, 子ノードを作る. 
