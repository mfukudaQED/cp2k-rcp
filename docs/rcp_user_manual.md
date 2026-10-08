# CP2K RCP 利用マニュアル（GPW・共線スピン）

**対象:** CP2K 2026.2 + RCP patch、開発ブランチ `rcp-development`

**更新:** 2026-10-08

> 本書は実際の計算・出力・検証のためのマニュアル。数式の詳細や実装履歴は [RCP詳細チュートリアル](rcp_tutorial.md) を参照。
> **非共線スピン・SOC・ZORA/DKHなどの相対論計算は保留中**であり、本書の使用対象外。
>
> **公開配布版について:** このリポジトリはCP2K本体から独立したpatch/例題配布パッケージ。
> `/path/to/cp2k-rcp` はこのリポジトリのclone先、`/path/to/patched-cp2k` は
> パッチを適用・ビルドした別のCP2Kソースツリーを意味する。
> 詳しい適用方法は[トップページ](../README.md)を参照。
> 大規模な5chain/diamond検証データや計算済みCubeは配布していない。

### このマニュアルの使い方

**すぐ計算したい場合は2・3節**（実行・入力）、設定の意味を調べる場合は**4・5節**、
Cubeの物理的意味と解析は**6～8節**、既存サンプルは**9節**、
問題が起きた場合は**10節**を参照する。

## 1. 対応機能と前提条件

| 計算条件 | 対応状況 | 補足 |
|---|---|---|
| QUICKSTEP / GPW | 対応 | 主にGTH擬ポテンシャル・PBEで検証 |
| Γ点・非周期分子 | 対応 | H₂、ベンゼンなど |
| Γ点・周期系 | 対応 | 5chainなど |
| 複数k点 | 対応 | 複素Bloch軌道、対称性縮約、MPIグループ |
| RKS（スピン非分極） | 対応 | スピン縮退2を考慮 |
| collinear UKS（スピン分極） | 対応 | α/βを合計して出力 |
| Fermi–Dirac smearing | 対応 | 300 Kの表面モデルで動作を確認。金属一般での精度保証ではない |
| GAPW / 全電子 / ZORA / DKH | **未対応** | RCPを要求すると明示的に停止 |
| 非共線・SOC・スピン螺旋 | **未対応** | 今回のマニュアル対象外 |
| Hybrid / HF | 未検証 | 別途検証が必要 |

RCPは**収束後のDFT電子状態を解析する機能**。SCF収束や構造最適化の手法を置き換えるものではない。

### 利用するバイナリ

パッチ適用・ビルド後の環境（パスは各自のCP2Kインストール先に変更）:

~~~bash
source /path/to/patched-cp2k/install/cp2k_env
which cp2k.psmp
~~~

標準配布のCP2K 2026.2に開発版RCP入力があるとは限らない。
**パッチ適用済みバイナリ**を使い、ビルド手順は[公開README](../README.md)と
CP2Kの公式ビルドガイドに従う。

## 2. shamで最初のRCP計算を実行する

以下はbash用コマンド。元のテストを上書きせず別の作業場所で実行する。

### 2.1 Γ点・H₂（RKS）

~~~bash
cd /path/to/cp2k-rcp
mkdir -p scratch/h2_gamma
cp examples/h2/H2-rcp.inp scratch/h2_gamma/
cd scratch/h2_gamma

sbatch --export=ALL,CP2K_ENV=/path/to/patched-cp2k/install/cp2k_env \
  ../../tools/run_cp2k_sham.slurm H2-rcp.inp h2_gamma.out
squeue -u "$USER"

# 終了後に確認
grep 'SCF run converged' h2_gamma.out
grep 'PROGRAM ENDED AT' h2_gamma.out
grep '^ RCP|' h2_gamma.out
ls -lh *.cube
~~~

### 2.2 複数k点・H₂（RKS）

~~~bash
cd /path/to/cp2k-rcp
mkdir -p scratch/h2_kpoints
cp examples/h2/H2-rcp-kpoints.inp scratch/h2_kpoints/
cd scratch/h2_kpoints
sbatch --export=ALL,CP2K_ENV=/path/to/patched-cp2k/install/cp2k_env \
  ../../tools/run_cp2k_sham.slurm H2-rcp-kpoints.inp h2_kpoints.out
~~~

この入力はMonkhorst–Pack `2 1 1`、全k点格子、複素Bloch軌道を使う。どちらも正常終了した既存回帰テストの完全入力である。

**注:** 2つのH₂回帰テストはファイル数を抑えるため `PRINT_DENSITY_WINDOW F` としてある。
そのため既定の実行では窓密度Cubeを作らず、**3種類のCube**を出力する。窓電子密度を解析したい場合は
コピー先の入力で `PRINT_DENSITY_WINDOW T` に変更する。`PRINT_RCP T` なら窓密度の数値積分は引き続き標準出力に表示される。

### 2.3 sham固有のMPI設定

`tools/run_cp2k_sham.slurm` はSlurm用バッチスクリプトで、MPI起動に `srun` ではなく `mpiexec` を使う。実際の実行部分:

~~~bash
export I_MPI_COLL_EXTERNAL=no
export I_MPI_FABRICS=shm:ofi
export OMP_NUM_THREADS=1
mpiexec -n "$SLURM_NTASKS" cp2k.psmp -i input.inp -o output.out
~~~

この起動方法は**Slurmジョブ内部で使用**する。ランク数・スレッド数・実行時間は通常のSlurm設定に従って調整する。上記はMPI 2ランク・1スレッドの単純な場合を想定した抜粋。

## 3. 入力ファイルの書き方

### 3.1 RCPを有効にする

通常のCP2K入力で `&FORCE_EVAL / &DFT / &PRINT` の中に `&RCP ON` を入れる。

~~~text
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    ! BASIS_SET_FILE_NAME、POTENTIAL_FILE_NAME、
    ! MGRID、SCF、XCなど通常のDFT設定をここに書く
    &PRINT
      &RCP ON
        ENERGY_LOWER [eV] -3.0
        ENERGY_UPPER [eV] 0.0
        BROADENING_LOWER [eV] 0.05
        BROADENING_UPPER [eV] 0.05
        DENSITY_CUTOFF 1.0E-7

        PRINT_DENSITY_WINDOW T
        PRINT_KINETIC_ENERGY_DENSITY T
        PRINT_REGIONAL_ENERGY_DENSITY T
        PRINT_RCP T

        MPI_IO F
        STRIDE 1 1 1
      &END RCP
    &END PRINT
  &END DFT
&END FORCE_EVAL
~~~

これは**RCPの挿入箇所を示す抜粋**であり単独では実行できない。完全入力は `examples/h2/`、`examples/benzene/`、`examples/c2h5/` にある。

### 3.2 RKS / UKS

共線UKSの例（doublet）:

~~~text
&DFT
  UKS T
  MULTIPLICITY 2
  ...
&END DFT
~~~

検証済みの開殻分子の完全入力: `examples/c2h5/c2h5.inp`。k点でUKS経路を通す回帰試験: `examples/h2/H2-rcp-kpoints-uks.inp`（UKS singlet）。

**RCPの4種類のCubeは、いずれもα/βスピン成分を合算した量。** スピン分解RCP Cubeの出力機能ではない。

### 3.3 Γ点 / 複数k点

Γ点専用であれば `&KPOINTS` を省略する。例えば4×4×1を使う場合:

~~~text
&KPOINTS
  SCHEME MONKHORST-PACK 4 4 1
  FULL_GRID OFF
  SYMMETRY ON
  WAVEFUNCTIONS COMPLEX
&END KPOINTS
~~~

- `FULL_GRID OFF` + `SYMMETRY ON`: 対称性を使ってk点を縮約する。
- `FULL_GRID ON` + `SYMMETRY OFF`: 参照計算として全格子を扱う。
- `PARALLEL_GROUP_SIZE 1`: k点MPIグループの動作検証例がある。
- 明示した `1 1 1` k点と、`&KPOINTS`なしのΓ点では内部処理が異なる。両経路のRCP整合性は数値的に検証済みだが、ビット単位の一致を仮定しない。

### 3.4 セル・擬ポテンシャル

分子系では `&CELL / PERIODIC NONE` と `&POISSON / PERIODIC NONE` を整合させる。スラブ系では、十分な真空を持つ3D周期セルに `k_x k_y 1` を設定する例がある。実際のPoisson解法・真空厚・双極子補正は通常のCP2K計算として適切に選択すること。

GPW/GTH擬ポテンシャルで求めるRCPの電子密度は**価電子密度**であり、GAPWの全電子密度とは異なる。


## 4. 各キーワードと推奨設定

下記は **`src/input_cp2k_print_dft.F` の定義を照合したデフォルト値**。すべて `&DFT / &PRINT / &RCP` の内側で設定する。

| キーワード | デフォルト | 役割 |
|---|---:|---|
| `ENERGY_LOWER` | −3.0 eV | 基準化学ポテンシャルから測った窓の下端 |
| `ENERGY_UPPER` | 0.0 eV | 同じく窓の上端 |
| `BROADENING_LOWER` | 0.001 eV | 下端のFermi関数の広がり |
| `BROADENING_UPPER` | 0.001 eV | 上端のFermi関数の広がり |
| `EPS_FILTER` | 1.0×10⁻¹⁴ | Γ点での疎行列演算のフィルタ閾値 |
| `DENSITY_CUTOFF` | 1.0×10⁻¹² a.u. | RCPの分母に用いる窓電子密度の閾値 |
| `PRINT_DENSITY_WINDOW` | T | 窓電子密度Cube |
| `PRINT_KINETIC_ENERGY_DENSITY` | T | 全占有電子のLaplacian型運動エネルギー密度Cube |
| `PRINT_REGIONAL_ENERGY_DENSITY` | T | 窓選択後の領域エネルギー密度Cube |
| `PRINT_RCP` | T | RCP Cube |
| `MPI_IO` | T | Cube出力のMPI-I/O |
| `STRIDE` | 1 1 1 | CubeのX/Y/Z方向の間引き |

エネルギー窓はデフォルトに頼らず、特に **`[eV]` を明記**する。

### 用途別の実用設定

| 用途 | 窓の例 | broadeningの例 | density cutoff | 出力 |
|---|---|---|---|---|
| 絶縁体・分子の狭い窓 | −3～0 eV | 0.001 eV | 1×10⁻¹² | 4種類、`STRIDE 1 1 1` |
| 表面のRCP解析 | −3～0 eV | 0.05 eV | 1×10⁻⁷ | 電子密度・RCP・領域エネルギー、必要なら運動エネルギー |
| near-metallic探索 | −1～0 eV | 0.01 eV | 系に応じて調整 | まず窓密度とRCPを同時出力 |
| 出力容量の節約 | 目的に応じる | 同上 | 同上 | `STRIDE 2 2 2`、不要なCubeをF |

これらは**設定例であり普遍的な推奨値ではない**。窓幅・broadening・閾値によるRCP変化も収束検証の対象になる。

- `ENERGY_LOWER` は必ず `ENERGY_UPPER` より小さくする。
- `BROADENING_LOWER` と `BROADENING_UPPER` は正にする。
- `PRINT_RCP T` の場合、分母である窓密度と分子である領域エネルギー密度は、各CubeをFにしても内部で計算される。
- `DENSITY_CUTOFF` を大きくすると真空・節近傍で数値の不安定な領域をより強くゼロ化する。**Cubeの可視化時だけに作用する閾値ではない。**
- `EPS_FILTER` は主にΓ点の密度行列構築・直交化で用いる閾値。k点経路の窓占有を制御するパラメータではない。
- `MPI_IO F` はMPI並列計算を無効化する意味ではない。**Cubeの出力方法**を切り替える設定である。
- `STRIDE` はCube出力の空間サンプリングを間引く。SCFや内部のRCP場そのものを粗いグリッドで再計算するものではない。

### 設定する前に決めること

**「どの電子状態を見たいか」**と**「RCPをどこに表示したいか」**を分けて考える。例えばダングリングボンドの空間分布と表面への吸着位置を比較したい場合、まず窓電子密度を確認して目的の軌道が十分含まれているかを見て、次にRCPと運動エネルギー密度を可視化する。

## 5. RCPの定義と窓占有

### 5.1 基準準位

RCPでは窓を**絶対KSエネルギーではなく基準化学ポテンシャルに対する相対値**で指定する。

$$
E_{\mathrm{lower}} \le \varepsilon_{n\mathbf k\sigma}-\mu \le E_{\mathrm{upper}}.
$$

- **SCF smearingなし:** 全スピン・全k点のHOMOとLUMOから、`μ=(ε_HOMO+ε_LUMO)/2` を求める。
- **SCF Fermi–Dirac smearingあり:** CP2KのSCFが求めたFermiエネルギーを `μ` とする。

したがって `ENERGY_LOWER [eV] -3.0` と `ENERGY_UPPER [eV] 0.0` は、基準準位からおよそ−3 eV～0 eVの軌道成分を選ぶことを意味する。

### 5.2 滑らかなエネルギー窓

各KS固有状態に対して、窓端を滑らかにしたFermi関数の差

$$
F(x;E,\delta)=\frac{1}{1+\exp[(x-E)/\delta]},\qquad
W_{n\mathbf k\sigma}
=\operatorname{clip}_{[0,1]}\left[
F(\varepsilon_{n\mathbf k\sigma}-\mu;E_{\mathrm{upper}},\delta_{\mathrm{upper}})
-F(\varepsilon_{n\mathbf k\sigma}-\mu;E_{\mathrm{lower}},\delta_{\mathrm{lower}})
\right]
$$

を窓占有の重みとして使う。`clip` は数値を0～1に収める処理である。

SCFの電子温度と、`BROADENING_LOWER/UPPER` **は別のもの**。窓の端に近い準位は部分的な窓占有となるため、窓電子数が整数にならなくても異常ではない。

### 5.3 出力する物理量

スピン合計の窓電子密度は

$$
n_{\mathrm{EW}}(\mathbf r)=
\sum_{\sigma,n,\mathbf k}
w_{\mathbf k}\,g_{\sigma}W_{n\mathbf k\sigma}\,
|\psi_{n\mathbf k\sigma}(\mathbf r)|^2
$$

である。RKSは `g=2`、UKSは各スピンチャネルで `g=1`。

窓密度行列 `P_EW` に対して、Laplacian型とgradient型の運動エネルギー密度を

$$
t_L^{\rm EW}(\mathbf r)
=-\frac14\sum_{\mu\nu}P^{\rm EW}_{\mu\nu}
\bigl[\phi_\mu\nabla^2\phi_\nu+(\nabla^2\phi_\mu)\phi_\nu\bigr],
$$

$$
t_G^{\rm EW}(\mathbf r)
=\frac12\sum_{\mu\nu}P^{\rm EW}_{\mu\nu}
\,\nabla\phi_\mu\cdot\nabla\phi_\nu
$$

と書くと、実装した領域エネルギー密度とRCPは

$$
\varepsilon_{\tau,\rm EW}(\mathbf r)
=-\frac12\left[t_L^{\rm EW}(\mathbf r)+t_G^{\rm EW}(\mathbf r)\right],
\qquad
\mu_R^\tau(\mathbf r)
=\frac{\varepsilon_{\tau,\rm EW}(\mathbf r)}{n_{\rm EW}(\mathbf r)}.
$$

**`n_EW ≤ DENSITY_CUTOFF` の格子点はRCPを0に設定する。**

`PRINT_KINETIC_ENERGY_DENSITY` が出力する `T_e` は、ここで使う窓選択後の `t_L^EW` ではなく、**全占有状態**に対応したLaplacian型運動エネルギー密度。RCPの分子と混同しないこと。


## 6. 計算結果・Cubeの読み方

`&GLOBAL / PROJECT mycalc` を指定し、RCPを出力した場合、代表的なファイル名は次のとおり。末尾の `1_0` などはCP2Kの反復番号に依存するため、固定文字列として解析スクリプトに埋め込まないこと。

| 出力キーワード | 代表的なCubeファイル | 物理量 | Cube値の単位 |
|---|---|---|---|
| `PRINT_DENSITY_WINDOW T` | `mycalc-RCP-DENSITY-WINDOW1_0.cube` | 窓電子密度 `n_EW(r)` | bohr⁻³ |
| `PRINT_KINETIC_ENERGY_DENSITY T` | `mycalc-RCP-KINETIC-ENERGY-DENSITY1_0.cube` | 全占有のLaplacian型 `T_e(r)` | Ha bohr⁻³ |
| `PRINT_REGIONAL_ENERGY_DENSITY T` | `mycalc-RCP-REGIONAL-ENERGY-DENSITY1_0.cube` | `ε_{τ,EW}(r)` | Ha bohr⁻³ |
| `PRINT_RCP T` | `mycalc-RCP1_0.cube` | `μ_R^τ(r)` | Ha |

RCPの単位は **Hartree**。eVに直す場合はCubeの数値に約 **27.211386** を掛ける。Cubeヘッダの格子ベクトルと原子座標の基本単位はbohrである（CP2Kの出力形式に依存するヘッダ仕様を確認すること）。

4種類を出す際は、同じSCF計算で作られるため共通の実空間グリッドを持つ。窓密度・領域エネルギー密度・RCPを組み合わせる際は、**同じ計算・同じ`STRIDE`・同じ反復番号**のファイルを使う。

### 6.1 典型的な使い方

- **RCPの空間分布を可視化:** `*-RCP1_0.cube` を可視化ソフトに読み込む。
- **電子的表面上でRCPを着色:** `T_e(r)` の等値面を作り、`μ_R^τ(r)` をカラーマップとして表示する。使用する等値面と色範囲は明記する。
- **比の信頼性を評価:** `n_EW(r)` も出力し、真空・軌道節など窓電子密度が極小の場所を解析から除く。
- **局所RCPと電子密度を比較:** Cubeを共通座標上で表示する。最大・最小値だけの比較は避ける。

### 6.2 Cube間の空白（viewer互換性）

現在のRCPブランチではCube数値の間に空白を入れる修正を適用済み。過去に作成したCubeが一部のviewerで読めない場合は、**再計算せずに**次のツールで空白のみを修正できる。

~~~bash
cd /path/to/cp2k-rcp
python3 tools/normalize_cube_spacing.py --check /path/to/cube_file.cube
python3 tools/normalize_cube_spacing.py --dry-run /path/to/cube_directory
python3 tools/normalize_cube_spacing.py /path/to/cube_directory
~~~

ディレクトリを指定した場合はその下のCubeを再帰的に処理する。既存データを扱うときは元のCubeのバックアップを確認する。

## 7. 標準出力の検証

通常、実行後に次を確認する。

~~~bash
grep 'SCF run converged' output.out
grep 'PROGRAM ENDED AT' output.out
grep '^ RCP|' output.out
~~~

**SCF未収束のRCP分布を物理的な結論に使わない。** 次の診断値が表示される（それぞれ対応する計算を有効にした場合）。

| 標準出力の項目（`RCP|`以下） | 解釈 |
|---|---|
| `Global HOMO`、`Global LUMO` | 窓の準位基準に利用するバンド端 |
| `HOMO-LUMO midpoint` | smearingなしでの基準 |
| `Reference chemical potential` | 実際に窓を置いた基準エネルギー（Ha） |
| `Relative energy window` | 入力で指定した相対窓。標準出力ではHa |
| `Electron count from sum(weights)` | KS状態の窓占有数とk点重みの総和 |
| `Electron count from Tr[P_window*S]` | 窓密度行列と重なり行列から求めた電子数 |
| `Electron count from grid integration` | 窓電子密度の実空間積分 |
| `\|Tr[P*S]-grid integral\|` | 実空間への変換・collocation・積分の整合性 |
| `Integral of Laplacian-form T_e` | 全占有Laplacian型運動エネルギーの実空間積分（Ha） |
| `CP2K kinetic energy` | CP2K内部での運動エネルギー（Ha） |
| `\|Integral(T_e)-E_kin\|` | 運動エネルギーの独立な整合性チェック |
| `Window Laplacian component integral` | 窓のLaplacian型運動エネルギー積分 |
| `Window gradient component integral` | 窓のgradient型運動エネルギー積分 |
| `Regional energy density integral` | 窓の領域エネルギー密度の体積積分 |

### 7.1 窓電子数の三重チェック

$$
N_{\rm EW}^{\rm weights}
\stackrel{?}{=}
\operatorname{Tr}[P_{\rm EW}S]
\stackrel{?}{=}
\int n_{\rm EW}(\mathbf r)\,d^3r .
$$

**3つの値の相互一致を確認**する。窓電子数が総価電子数と一致する必要はない。

検証済みH₂の例では、窓をほぼ全価電子状態を含むように広く取ったとき、

~~~text
RCP| Electron count from grid integration: 2.000000002782
RCP| Regional energy density integral [a.u.]: -1.134090736496
~~~

となった。ただし**数値は指定した窓、構造、基底、k点、計算条件に依存**する。別の系の値が同じになることを要求しない。

### 7.2 運動エネルギーの二重チェック

$$
\int T_e(\mathbf r)\,d^3r \stackrel{?}{=} E_{\rm kinetic}^{\rm CP2K}.
$$

また完全な積分では、適切な境界条件のもとで窓Laplacian成分と窓gradient成分が一致するはず。異常に大きな差が出る場合は、SCF、グリッド、境界条件、k点変換を調べる。

### 7.3 診断値が出ない場合

`&RCP ON` の指定、`&DFT / &PRINT` の階層、`PRINT_LEVEL`、SCF正常終了、利用中のバイナリのバージョンを確認する。`RCP|`の行があること自体はSCF収束の代用ではない。

## 8. k点収束と表面系の解析

### 8.1 k点数を変えるときに固定する条件

- **同じ原子配置とセル**（途中で再最適化しない）。
- 同じ汎関数・基底・擬ポテンシャル・`CUTOFF`・`REL_CUTOFF`。
- 同じ`SCF EPS_SCF`、smearing、窓上下端、窓のbroadening。
- 同じ実空間グリッド、`STRIDE`、`DENSITY_CUTOFF`。
- 同じ観察面と、同じ窓電子密度マスク。

比較する指標は**全エネルギーだけではなく**、窓電子数、窓密度の分布、RCPの断面・等値面・RMS差など。異なるk点数の計算間では同じ座標系での比較が前提。

### 8.2 真空・節領域を除く定量比較

窓密度 `n_EW` の極小領域で `μ_R^τ=ε_{τ,EW}/n_EW` は不安定になりやすい。基準計算の窓密度からマスク `Ω` を作り、**全メッシュで同じ格子点集合**を使って比較する。

$$
\Delta_{\rm RMS}(\mathcal K,\mathcal K_{\rm ref})
=
\sqrt{\frac{1}{|\Omega|}
\sum_{\mathbf r_i\in\Omega}
\left[
\mu_R^{\tau,\mathcal K}(\mathbf r_i)
-\mu_R^{\tau,\mathcal K_{\rm ref}}(\mathbf r_i)
\right]^2}.
$$

マスク閾値（例えば最大窓密度の0.1%、5%）を変えた結果も確認する。最も細かいk点メッシュは**暫定的な参照**であり、そのメッシュ自身とのRMS差0は収束の証明ではない。

### 8.3 ダイヤモンド(001)の既存検証

開発環境では表面スラブをΓ点から8×8のk点メッシュで評価した。
ただし、その計算入力・Cube・解析グラフと元の構造データは**この公開パッケージに含まれない**。
以下の数値は方法論上の参考であり、このパッケージ単体からの再現を保証しない。

この計算では6×6と8×8の**全エネルギー差は約1.84 meV/セル**まで低下したが、表面RCPのRMS差は**約3.89 eV**残っていた。**RCPのk点収束は8×8まででは確定していない。** このモデルはC₁₀H₄の固定幾何で、底面のH終端が近似的なので、実物質の定量予測にはスラブ厚や緩和構造も再検討すること。


## 9. 既存サンプルと検証データ

公開リポジトリルート `/path/to/cp2k-rcp` からの相対パス。

| 用途 | 入力または検証データ |
|---|---|
| H₂、非周期Γ点RKS | `examples/h2/H2-rcp.inp` |
| H₂、周期Γ点RKS | `examples/h2/H2-rcp-gamma-periodic.inp` |
| H₂、明示的1×1×1 k点 | `examples/h2/H2-rcp-kpoints-1x1x1.inp` |
| H₂、2×1×1 k点（全格子） | `examples/h2/H2-rcp-kpoints.inp` |
| H₂、2×1×1 k点（対称性縮約） | `examples/h2/H2-rcp-kpoints-sym.inp` |
| H₂、UKS singletで複数k点 | `examples/h2/H2-rcp-kpoints-uks.inp` |
| H₂、k点MPI並列グループ | `examples/h2/H2-rcp-kpoints-parallel.inp` |
| ベンゼン分子、RKS | `examples/benzene/benzene.inp` |
| C₂H₅ラジカル、UKS doublet | `examples/c2h5/c2h5.inp` |

### 9.1 回帰試験の実行結果を確認する

H₂の各k点入力を計算して、必要な`*.out`と`*.cube`を1つの
`regression-results/` ディレクトリに集めた場合に、SCF診断値とRCP Cubeを再検証できる。
**検証結果ファイル自体は配布していない。**

~~~bash
cd /path/to/cp2k-rcp
python3 tools/check_kpoint_consistency.py regression-results
~~~

既存の数値基準は `examples/h2/TEST_FILES.toml` に記録されている。

### 9.2 自分の計算を保存するとき

計算条件が分からなくならないよう、少なくとも次を一緒に保存する。

- CP2Kの入力 `*.inp`、標準出力 `*.out`、SCF収束情報。
- CP2Kのバイナリ/commitとパッチの版。
- 基底・擬ポテンシャル、汎関数、格子・原子配置。
- k点メッシュ、Fermi smearing、窓上下端・窓のbroadening。
- Cube、`DENSITY_CUTOFF`、`STRIDE`、可視化・比較で使ったマスク条件。

## 10. トラブルシューティング

| 症状 | 確認と対処 |
|---|---|
| `RCP` が入力として認識されない | パッチを適用した `cp2k.psmp` を使用しているか確認。通常版CP2Kと混同しない。 |
| `SCF run NOT converged` | まずSCFを収束させる。`EPS_SCF`、`MAX_SCF`、`MIXING`、`ADDED_MOS`、smearingなどを通常のDFT設定として見直す。 |
| `RCP ENERGY_LOWER must be smaller ...` | `ENERGY_LOWER < ENERGY_UPPER` を指定する。 |
| `RCP broadenings must be positive` | `BROADENING_LOWER / UPPER` を正にする。 |
| `RCP ... needs occupied and unoccupied bands` | k点経路で占有/非占有バンドが必要。`ADDED_MOS` や固有状態の本数を確認する。 |
| RCP Cubeの真空中に巨大値が出る | 窓電子密度と`DENSITY_CUTOFF`を確認。低密度領域を定量比較から除外する。 |
| Cubeが大きすぎる | `STRIDE 2 2 2`、不要な`PRINT_*`を`F`、分析対象領域を絞る。 |
| Cube I/Oで進まない、viewerで読めない | `MPI_IO F` を試す。過去のファイルは `normalize_cube_spacing.py` で確認。 |
| 窓電子数の3通りの値が合わない | SCF、同一k点設定、数値グリッド、積分・行列変換を確認。 |
| `GAPW ... one-center ... not implemented` | 現行RCPはGPW対象。GAPWの全電子RCPを意味する結果は生成できない。 |
| `RCP ... relativistic kinetic operator ... not implemented` | 明示的な相対論的運動エネルギー演算子は未対応。 |
| k点結果がΓ点と大きく違う | 原子配置・周期境界・k点格子の違いを確認し、窓電子数とRCP分布のk点収束を調べる。 |
| `srun` でMPIエラー | shamでは`mpiexec`を使用し、指定の`I_MPI_*`環境変数を設定する。 |

なお、`DENSITY_CUTOFF` を調節するだけで大きな差を隠すのは避ける。窓密度と領域エネルギー密度の両方を確認し、収束が悪い原因を切り分ける。

## 11. 実行前チェックリスト

1. [ ] CP2K 2026.2の**RCP実装済み**バイナリを使う。
2. [ ] GPW、RKSまたは共線UKSの入力になっている。
3. [ ] セル・周期性・Poisson解法・基底・擬ポテンシャルが通常のDFTとして適切。
4. [ ] Γ点か複数k点かを明示し、k点収束を検討した。
5. [ ] SCFが収束し、必要な非占有バンド数が確保されている。
6. [ ] 窓の基準 `μ`、`ENERGY_LOWER/UPPER`、broadeningの意味を把握した。
7. [ ] RCPと同時に`PRINT_DENSITY_WINDOW T`を使い、マスクの妥当性を確認する。
8. [ ] `PRINT_KINETIC_ENERGY_DENSITY T`が出すのは**全占有**の `T_e` と理解した。
9. [ ] `RCP|`の電子数・運動エネルギーの診断値をチェックした。
10. [ ] RCPの単位（Ha）、閾値、Cubeの`STRIDE`、可視化面・色範囲を記録した。
11. [ ] k点数だけでなくRCP分布そのものの収束を確認した。
12. [ ] shamでMPIを使う場合は`mpiexec`と必要な環境変数を使う。

---

**さらに詳しい資料**

- [RCP詳細チュートリアル](rcp_tutorial.md) — 定義、各種分子・周期系検証、Cubeの可視化
- `examples/h2/` — 自動回帰試験入力・基準値
- ダイヤモンド(001)の検証データは別途管理（公開パッケージに未同梱）
- [公開README](../README.md) — 導入・パッチ適用手順
- パッチ適用後のCP2K `src/input_cp2k_print_dft.F` — `&PRINT / &RCP`のキーワード定義
- パッチ適用後のCP2K `src/qs_energy_window.F` — RCP本体およびΓ点/k点密度行列経路
