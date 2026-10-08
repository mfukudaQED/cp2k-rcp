# CP2K RCP Tutorial

> **公開配布版の注記:** 本書は開発環境での検証記録も含む詳細チュートリアル。
> `DEV_TREE/` と表記する開発ソース・`local-sham` 検証データはこの公開リポジトリに含まれない。
> 同梱の入力は[サンプル一覧](../examples/README.md)、実行可能な公開手順は
> [README](../README.md) と [利用者マニュアル](rcp_user_manual.md) を参照。


> **利用者向け入口:** [CP2K RCP 利用マニュアル](rcp_user_manual.md) を先に参照。
> この文書は詳しい背景説明・実装経緯・個別検証例を残した詳細版。


## 1. 概要

この文書では、CP2K の rcp-development branch に実装した Regional Chemical Potential (RCP) を計算する方法を説明する。

現在、実装・検証済みの主な対象は次の通り。

- GPW
- Γ点計算、一般 k-point sampling（GPW）
- collinear RKS / UKS
- 分子系
- 周期系
- 整数占有
- Fermi–Dirac smearing を用いた near-metallic 系

RCP は

$$
\mu_R^\tau(\mathbf r)
=
\frac{\varepsilon_{\tau,\mathrm{ew}}(\mathbf r)}
     {n_{\mathrm{ew}}(\mathbf r)}
$$

で定義する。

ここで、

- \(n_{\mathrm{ew}}(\mathbf r)\): energy window に含まれる電子密度
- \(\varepsilon_{\tau,\mathrm{ew}}(\mathbf r)\): energy-window regional energy density
- \(\mu_R^\tau(\mathbf r)\): regional chemical potential

である。

この実装では、全占有状態に対する Laplacian-form kinetic-energy density

$$
T_e(\mathbf r)
=
-\frac14
\sum_{\mu\nu}P_{\mu\nu}
\left[
\phi_\mu(\mathbf r)\nabla^2\phi_\nu(\mathbf r)
+
\nabla^2\phi_\mu(\mathbf r)\phi_\nu(\mathbf r)
\right]
$$

も同時に計算できる。

---

## 2. RCP を有効にする最小入力

RCP の設定は

DFT → PRINT → RCP

の階層に書く。

典型例:

~~~text
&DFT
  ...
  &PRINT
    &RCP ON
      ENERGY_LOWER [eV] -3.0
      ENERGY_UPPER [eV]  0.0

      BROADENING_LOWER [eV] 0.001
      BROADENING_UPPER [eV] 0.001

      DENSITY_CUTOFF 1.0E-12

      PRINT_DENSITY_WINDOW T
      PRINT_KINETIC_ENERGY_DENSITY T
      PRINT_REGIONAL_ENERGY_DENSITY T
      PRINT_RCP T

      MPI_IO T
      STRIDE 1 1 1
    &END RCP
  &END PRINT
&END DFT
~~~

energy window は、絶対 KS energy ではなく **基準 chemical potential に対する相対値**で指定する。

例えば

~~~text
ENERGY_LOWER [eV] -3.0
ENERGY_UPPER [eV]  0.0
~~~

は、おおよそ

$$
\mu-3\ {\rm eV}
<
\varepsilon
<
\mu
$$

の範囲を選択する。

---

### 2.1 Γ点以外の k-point sampling

CP2K 2026.2 RCP patch の k-point 版では、SCF と同じ `&KPOINTS` を指定する。
例えば、2×1×1 Monkhorst–Pack sampling は次のように指定する。

~~~text
&DFT
  &KPOINTS
    SCHEME MONKHORST-PACK 2 1 1
    FULL_GRID ON
    SYMMETRY OFF
    WAVEFUNCTIONS COMPLEX
  &END KPOINTS
  &PRINT
    &RCP ON
      ENERGY_LOWER [eV] -20.0
      ENERGY_UPPER [eV] 1.0
      PRINT_DENSITY_WINDOW T
      PRINT_KINETIC_ENERGY_DENSITY T
      PRINT_REGIONAL_ENERGY_DENSITY T
      PRINT_RCP T
    &END RCP
  &END PRINT
&END DFT
~~~

k-point 計算は、各 \(\mathbf k\) の複素 Bloch 軌道と固有値を使う。
エネルギー窓密度は、CP2K と同じ正規化済み重み \(w_{\mathbf k}\) により

$$
n_{\mathrm{ew}}(\mathbf r)=
\sum_{\sigma}\sum_{\mathbf k} w_{\mathbf k}
\sum_n f_{n\mathbf k\sigma}^{\mathrm{ew}}
\left|\psi_{n\mathbf k\sigma}(\mathbf r)\right|^2
$$

を評価する。スピン縮退した RKS の係数は2、UKSではスピンごとに1。
固有状態の一時的な窓占有数から密度行列を構成し、標準の
`kpoint_density_transform` により Bloch 位相、結晶対称性と実空間 lattice images を扱う。
標準 SCF の占有数・k点密度行列は RCP 出力の前後で復元される。

テスト入力例:
`tests/QS/regtest-rcp/H2-rcp-kpoints.inp` および
`tests/QS/regtest-rcp/H2-rcp-kpoints-sym.inp`。
精度検証は通常と同様に、\(\mathrm{Tr}[P_{\mathrm{ew}}S]\)、
密度の実空間積分、kinetic energy の積分値を比較する。

---

## 3. RCP キーワード

| keyword | default | 説明 |
|---|---:|---|
| ENERGY_LOWER | -3.0 eV | energy window 下端。chemical potential からの相対値 |
| ENERGY_UPPER | 0.0 eV | energy window 上端 |
| BROADENING_LOWER | 0.001 eV | window 下端の Fermi broadening |
| BROADENING_UPPER | 0.001 eV | window 上端の Fermi broadening |
| EPS_FILTER | 1.0E-14 | sparse matrix operation の filtering threshold |
| DENSITY_CUTOFF | 1.0E-12 | RCP の分母 \(n_{\mathrm{ew}}\) に対する cutoff |
| PRINT_DENSITY_WINDOW | T | \(n_{\mathrm{ew}}\) cube を出力 |
| PRINT_KINETIC_ENERGY_DENSITY | T | \(T_e\) cube を出力 |
| PRINT_REGIONAL_ENERGY_DENSITY | T | \(\varepsilon_{\tau,\mathrm{ew}}\) cube を出力 |
| PRINT_RCP | T | \(\mu_R^\tau\) cube を出力 |
| MPI_IO | T | cube output に parallel MPI-I/O を使う |
| STRIDE | 1 1 1 | cube grid の X/Y/Z stride |

エネルギー関係の値は単位を明示することを推奨する。

~~~text
ENERGY_LOWER [eV] -3.0
ENERGY_UPPER [eV] 0.0
BROADENING_LOWER [eV] 0.001
BROADENING_UPPER [eV] 0.001
~~~

---

## 4. 基準 chemical potential

### 4.1 smearing を使わない場合

通常の整数占有計算では、

$$
\mu
=
\frac{
\varepsilon_{\rm HOMO}
+
\varepsilon_{\rm LUMO}
}{2}
$$

を RCP energy window の基準とする。

UKS では、全 spin channel から

- global HOMO = 各 spin channel の HOMO の最大値
- global LUMO = 各 spin channel の LUMO の最小値

を求め、その midpoint を用いる。

標準出力には次の情報が出る。

~~~text
RCP| Global HOMO [a.u.]
RCP| Global LUMO [a.u.]
RCP| HOMO-LUMO midpoint [a.u.]
RCP| Reference chemical potential [a.u.]
~~~

smearing 無しでは、

Reference chemical potential = HOMO-LUMO midpoint

となる。

### 4.2 Fermi–Dirac smearing を使う場合

SCF で smearing が有効な場合は、midpoint ではなく SCF が求めた Fermi energy を RCP の基準 chemical potential に使う。

~~~text
&SCF
  ...
  &SMEAR
    METHOD FERMI_DIRAC
    ELECTRONIC_TEMPERATURE [K] 300
  &END SMEAR
&END SCF
~~~

この場合、

Reference chemical potential = SCF Fermi energy

となる。

metallic / near-metallic 系ではこちらを推奨する。

---

## 5. Energy window の Fermi weight

energy window は完全な step function ではなく、Fermi function の差で滑らかに定義される。

概念的には各 KS state に対して

$$
w_i
=
f(
\varepsilon_i-\mu;
E_{\rm upper},
\sigma_{\rm upper}
)
-
f(
\varepsilon_i-\mu;
E_{\rm lower},
\sigma_{\rm lower}
)
$$

を使う。

したがって、window edge にちょうど位置する state は broadening が有限なら部分占有になる。

例えば finite molecule の単体テストで HOMO をほぼ完全に含めたい場合には、ENERGY_UPPER を 0 eV より少し正側にする方法もある。

一方、表面・金属系では

~~~text
ENERGY_UPPER [eV] 0.0
~~~

として Fermi level 近傍を滑らかに選ぶ使い方が自然である。

---

## 6. Closed-shell 分子の例: benzene

検証済みの完全入力:

~~~text
DEV_TREE/local-sham/rcp_validation/benzene/input/benzene.inp
~~~

RCP 部分:

~~~text
&PRINT
  &RCP ON
    ENERGY_LOWER [eV] -3.0
    ENERGY_UPPER [eV] 0.0

    BROADENING_LOWER [eV] 0.001
    BROADENING_UPPER [eV] 0.001

    DENSITY_CUTOFF 1.0E-12

    PRINT_DENSITY_WINDOW T
    PRINT_KINETIC_ENERGY_DENSITY T
    PRINT_REGIONAL_ENERGY_DENSITY T
    PRINT_RCP T

    STRIDE 1 1 1
  &END RCP
&END PRINT
~~~

非周期分子では CELL と POISSON の周期性を合わせる。

~~~text
&POISSON
  PERIODIC NONE
  PSOLVER WAVELET
&END POISSON

...

&CELL
  ...
  PERIODIC NONE
&END CELL
~~~

---

## 7. Open-shell UKS の例: C2H5

現在の RCP は collinear UKS に対応している。

検証済みの完全入力:

~~~text
DEV_TREE/local-sham/rcp_validation/c2h5/input/c2h5.inp
~~~

C2H5 radical では

~~~text
&DFT
  UKS T
  MULTIPLICITY 2
  ...
~~~

とする。

RCP block:

~~~text
&PRINT
  &RCP ON
    ENERGY_LOWER [eV] -2.0
    ENERGY_UPPER [eV] 0.0

    BROADENING_LOWER [eV] 0.001
    BROADENING_UPPER [eV] 0.001

    DENSITY_CUTOFF 1.0E-12

    PRINT_DENSITY_WINDOW T
    PRINT_KINETIC_ENERGY_DENSITY T
    PRINT_REGIONAL_ENERGY_DENSITY T
    PRINT_RCP T

    STRIDE 1 1 1
  &END RCP
&END PRINT
~~~

RCP cube は alpha / beta を別々には出さず、全 collinear spin channel を合算した実空間場を出力する。

---

## 8. 周期系・near-metallic 系の例: 5chain

検証済みの完全入力:

~~~text
DEV_TREE/local-sham/rcp_validation/5chain/input/5chain_smear_serial_full.inp
~~~

5chain は frontier gap が約 0.025 eV で、300 K の \(k_BT\) とほぼ同程度だった。

そのため、整数占有 OT ではなく diagonalization + Fermi–Dirac smearing を用いた。

SCF の主要部分:

~~~text
&SCF
  EPS_SCF 1.0E-7
  MAX_SCF 100
  SCF_GUESS RESTART

  ADDED_MOS 50

  &MIXING
    ALPHA 0.20
    BETA 1.0
    METHOD BROYDEN_MIXING
    NBROYDEN 12
  &END MIXING

  &SMEAR
    METHOD FERMI_DIRAC
    ELECTRONIC_TEMPERATURE [K] 300
  &END SMEAR
&END SCF
~~~

RCP 部分:

~~~text
&PRINT
  &RCP ON
    ENERGY_LOWER [eV] -1.0
    ENERGY_UPPER [eV] 0.0

    BROADENING_LOWER [eV] 0.01
    BROADENING_UPPER [eV] 0.01

    DENSITY_CUTOFF 1.0E-12

    PRINT_DENSITY_WINDOW F
    PRINT_KINETIC_ENERGY_DENSITY T
    PRINT_REGIONAL_ENERGY_DENSITY F
    PRINT_RCP T

    MPI_IO F
    STRIDE 1 1 1
  &END RCP
&END PRINT
~~~

この例では大規模 periodic cube の parallel MPI-I/O が停止したため、最終計算では

~~~text
MPI_IO F
~~~

を使用した。

大規模系で cube output が進まない場合は、まず MPI_IO F を試すことを推奨する。

---

## 9. どの cube を出力するべきか

### 9.1 実装検証をしたい場合

4種類すべて出す。

~~~text
PRINT_DENSITY_WINDOW T
PRINT_KINETIC_ENERGY_DENSITY T
PRINT_REGIONAL_ENERGY_DENSITY T
PRINT_RCP T
~~~

### 9.2 表面の RCP 可視化だけが目的の場合

通常は \(T_e\) と RCP の2つで十分。

~~~text
PRINT_DENSITY_WINDOW F
PRINT_KINETIC_ENERGY_DENSITY T
PRINT_REGIONAL_ENERGY_DENSITY F
PRINT_RCP T
~~~

PRINT_RCP T の場合、RCP を作るために \(n_{\mathrm{ew}}\) と \(\varepsilon_{\tau,\mathrm{ew}}\) は内部では自動的に計算される。

したがって cube として保存しないなら、それぞれの PRINT flag を F にしてよい。

---

## 10. 出力 cube ファイル

GLOBAL section で

~~~text
PROJECT mycalc
~~~

とした場合、代表的なファイル名は次のようになる。

### 10.1 Energy-window electron density

~~~text
mycalc-RCP-DENSITY-WINDOW1_0.cube
~~~

内容:

$$
n_{\mathrm{ew}}(\mathbf r)
$$

単位:

$$
{\rm bohr}^{-3}
$$

RKS では closed-shell spin degeneracy 2 を含む電子数密度である。

この cube を全空間積分すると、window に含まれる電子数に対応する。

### 10.2 Laplacian-form kinetic-energy density

~~~text
mycalc-RCP-KINETIC-ENERGY-DENSITY1_0.cube
~~~

内容:

$$
T_e(\mathbf r)
=
-\frac14
\sum_{\mu\nu}P_{\mu\nu}
\left[
\phi_\mu\nabla^2\phi_\nu
+
(\nabla^2\phi_\mu)\phi_\nu
\right]
$$

単位:

$$
{\rm Hartree}\,{\rm bohr}^{-3}
$$

注意:

これは positive-definite kinetic-energy density \(\tau(\mathbf r)\) とは異なる。

### 10.3 Energy-window regional energy density

~~~text
mycalc-RCP-REGIONAL-ENERGY-DENSITY1_0.cube
~~~

内容:

$$
\varepsilon_{\tau,\mathrm{ew}}(\mathbf r)
=
\frac18
\sum_{\mu\nu}
P^{\mathrm{ew}}_{\mu\nu}
\left[
\phi_\mu\nabla^2\phi_\nu
+
(\nabla^2\phi_\mu)\phi_\nu
\right]
-
\frac14
\sum_{\mu\nu}
P^{\mathrm{ew}}_{\mu\nu}
\nabla\phi_\mu\cdot\nabla\phi_\nu
$$

単位:

$$
{\rm Hartree}\,{\rm bohr}^{-3}
$$

### 10.4 Regional chemical potential

~~~text
mycalc-RCP1_0.cube
~~~

内容:

$$
\mu_R^\tau(\mathbf r)
=
\frac{
\varepsilon_{\tau,\mathrm{ew}}(\mathbf r)
}{
n_{\mathrm{ew}}(\mathbf r)
}
$$

単位:

$$
{\rm Hartree}
$$

cube title にも

~~~text
RCP REGIONAL CHEMICAL POTENTIAL [HARTREE]
~~~

と出る。

ファイル名の 1_0 などの suffix は CP2K の iteration 番号に依存する。

スクリプトでは

~~~text
*-RCP*.cube
~~~

のように検索する方が安全である。

---

## 11. DENSITY_CUTOFF

RCP は

$$
\mu_R^\tau
=
\frac{
\varepsilon_{\tau,\mathrm{ew}}
}{
n_{\mathrm{ew}}
}
$$

という比なので、vacuum や orbital node では \(n_{\mathrm{ew}}\to 0\) となり、巨大値や数値ノイズが発生しやすい。

そのため、

~~~text
DENSITY_CUTOFF 1.0E-12
~~~

を設定すると、

$$
n_{\mathrm{ew}}(\mathbf r)
\le
10^{-12}
$$

の grid point の RCP を 0 にする。

これは表示上の clip ではなく、RCP ratio を計算するときの denominator mask である。

volume 全体の correlation、min/max、histogram などを計算するときも、低 density 領域を物理的な比較領域から除外することを推奨する。

---

## 12. STRIDE

~~~text
STRIDE 1 1 1
~~~

は CP2K real-space grid の全点を cube に書く。

例えば

~~~text
STRIDE 2 2 2
~~~

では X/Y/Z 各方向で2点ごとに出力するため、cube の grid point 数は概ね 1/8 になる。

推奨:

- 実装検証・定量比較: STRIDE 1 1 1
- 大規模系の可視化だけ: STRIDE 2 2 2 以上も検討

---

## 13. MPI_IO

default:

~~~text
MPI_IO T
~~~

小～中規模系では parallel cube output が利用できる。

一方、5chain の大規模 periodic test では parallel MPI-I/O が停止したため、

~~~text
MPI_IO F
~~~

で serial cube output を使った。

計算本体が終了しているのに cube file が 0 byte のまま増えない場合などは MPI_IO F を試す。

---

## 14. 標準出力で確認すべき診断値

RCP が実行されると、CP2K output に RCP| で始まる診断値が表示される。

代表例:

~~~text
RCP| Number of collinear spin channels:
RCP| Global HOMO [a.u.]:
RCP| Global LUMO [a.u.]:
RCP| HOMO-LUMO midpoint [a.u.]:
RCP| Reference chemical potential [a.u.]:

RCP| Relative energy window [a.u.]:
RCP| Lower/upper broadening [a.u.]:

RCP| Electron count from sum(weights):
RCP| Electron count from Tr[P_window*S]:
RCP| Electron count from grid integration:
RCP| |Tr[P*S]-grid integral|:

RCP| Integral of Laplacian-form T_e [a.u.]:
RCP| CP2K kinetic energy [a.u.]:
RCP| |Integral(T_e)-E_kin| [a.u.]:

RCP| Window Laplacian component integral [a.u.]:
RCP| Window gradient component integral [a.u.]:
RCP| Regional energy density integral [a.u.]:

RCP| Density cutoff for RCP ratio [a.u.]:
~~~

---

## 15. Window density の consistency check

次の3種類の electron count が一致することを確認する。

Fermi weights:

$$
N_{\mathrm{ew}}
=
\sum_i g_s w_i
$$

density matrix:

$$
N_{\mathrm{ew}}
=
{\rm Tr}
\left[
P^{\mathrm{ew}}S
\right]
$$

real-space density:

$$
N_{\mathrm{ew}}
=
\int
n_{\mathrm{ew}}(\mathbf r)
d^3r
$$

対応する出力:

~~~text
Electron count from sum(weights)
Electron count from Tr[P_window*S]
Electron count from grid integration
~~~

また、

~~~text
|Tr[P*S]-grid integral|
~~~

が十分小さいことを確認する。

---

## 16. Kinetic-energy density の consistency check

Laplacian-form kinetic-energy density は

$$
\int
T_e(\mathbf r)
d^3r
=
E_{\rm kin}
$$

を満たすべきである。

したがって、

~~~text
RCP| Integral of Laplacian-form T_e [a.u.]
RCP| CP2K kinetic energy [a.u.]
RCP| |Integral(T_e)-E_kin| [a.u.]
~~~

を確認する。

最後の差が十分小さければ、real-space derivative collocation の強い内部チェックになる。

---

## 17. 実行方法

### 17.1 一般的な実行

~~~bash
source DEV_TREE/install/cp2k_env

mpiexec -n 8 cp2k.psmp \
  -i input.inp \
  -o output.out
~~~

### 17.2 sham

sham 用 wrapper:

~~~text
DEV_TREE/local-sham/run_cp2k.slurm
~~~

例:

~~~bash
cd DEV_TREE

sbatch -n 8 \
  local-sham/run_cp2k.slurm \
  input.inp \
  output.out
~~~

sham では MPI 実行に srun ではなく mpiexec を使う。

---

## 18. RCP の可視化

表面系では、RCP の3D volume 全体を見るよりも、電子的 interface 上に RCP を表示する方が有用な場合が多い。

これまでの検証では、

$$
T_e(\mathbf r)
=
10^{-5}
$$

の isosurface を作り、その surface を

$$
\mu_R^\tau(\mathbf r)
$$

で着色した。

その場合に必要なファイルは

~~~text
*-RCP-KINETIC-ENERGY-DENSITY*.cube
*-RCP*.cube
~~~

である。

5chain の OpenMX 比較では例として

$$
-0.4
\le
\mu_R^\tau
\le
-0.25
\ {\rm Ha}
$$

を color range として使用した。

validation 済みの可視化用 cube は

~~~text
DEV_TREE/local-sham/rcp_validation/visualize/
~~~

にまとめてある。

代表例:

~~~text
benzene_RCP.cube
benzene_Te.cube

c2h5_RCP.cube
c2h5_Te.cube

5chain_RCP.cube
5chain_Te.cube
~~~

---

## 19. Cube file の whitespace compatibility

現在の branch では cube の volumetric values 間に明示的な whitespace を入れて出力する。

例えば、

~~~text
  0.00000E+000 -0.19538E+001 -0.20544E+001
~~~

のようになる。

この形式なら、固定幅 Fortran field を解釈しない viewer でも読みやすい。

古い cube の whitespace を修正する converter:

~~~text
DEV_TREE/local-sham/tools/normalize_cube_spacing.py
~~~

確認のみ:

~~~bash
python3 \
  DEV_TREE/local-sham/tools/normalize_cube_spacing.py \
  --check old.cube
~~~

変換:

~~~bash
python3 \
  DEV_TREE/local-sham/tools/normalize_cube_spacing.py \
  old.cube
~~~

ディレクトリを指定すると、その下の cube を再帰的に処理する。

converter は数値を再計算せず、元の数値 token を保持したまま whitespace だけを正規化する。

---

## 20. 現在の制限事項

### 20.1 k-point sampling

GPW の Γ点と一般 k-point sampling に対応する。
複素 Bloch 軌道および time-reversal/inversion symmetry を伴う
k点縮約は CP2K 標準の密度行列変換ルーチンを使用する。
現在の小規模回帰テストでは、1×1×1、2×1×1、対称性あり/なしを対象とする。
大規模系・金属系・異なる k-point mesh に対しては、独立した収束検証が必要。

### 20.2 collinear spin のみ

CP2K 2026.2 のRCP Cube出力が現在対応しているのは、RKSまたはcollinear UKSのみである。
SOCを使ったpost-SCFスピノル解析や、自己無撞着なnoncollinear磁性を
このRCP計算が自動的に扱うわけではない。

スピノルRCP計算核（電荷、磁化3成分、非相対論的なkinetic/regional energy）
は `src/rcp_spinor.F` にあり、単体テストは
`tests/QS/regtest-rcp/test_spinor_kernel.f90` にある。
ただしスピノル固有状態の読み取りと実空間collocationは未実装である。

GAPWでのZORA/DKHは、GPW用RCPの密度collocatorでは全電子補正を扱えない。
対応前に無効な結果を出すのを防ぐため、RCP計算を明示的に中断する。
対応方針は `docs/methods/rcp_spinor_relativistic_design.md` を参照する。

対応:

- RKS
- collinear UKS

現時点で未対応:

- noncollinear spin
- noncollinear spinor wavefunctions（複素 Bloch 軌道は対応）
- spin spiral

### 20.3 現在十分に検証している electronic-structure setting

検証済みの中心は

- GPW
- GTH pseudopotential
- semilocal PBE
- Γ point、および小規模 GPW k-point 系

である。

GAPWは現在RCP出力に対応せず、実行時に明示的に停止する。
hybrid/HF、一般的な大規模k点mesh、金属状態一般のk点smearingに対する系統的検証は未完了である。
一方、ダイヤモンド(001)スラブのPBE・300 K smearingでは、Γ点から8×8までの
k点RCPを実行済みである（`local-sham/diamond001_kconvergence/analysis/REPORT.md`）。
この検証では全エネルギーがほぼ収束してもRCP分布には差が残った。

### 20.4 大規模系

RCP は SCF 後に

- Γ点経路では overlap matrix の \(S^{-1/2}\) と dense KS diagonalization
- k点経路では SCF の Bloch MO を再利用して energy-window occupation を構成
- energy-window density matrix 構築
- real-space derivative collocation
- cube output

を行う。

したがって AO 数が大きい系では、通常の SCF より

- メモリ
- dense linear algebra
- I/O

が重くなる。

大規模系では

~~~text
MPI_IO F
~~~

および

~~~text
STRIDE 2 2 2
~~~

を必要に応じて検討する。

---

## 21. 推奨 RCP block

### 21.1 分子・絶縁体

~~~text
&RCP ON
  ENERGY_LOWER [eV] -3.0
  ENERGY_UPPER [eV] 0.0

  BROADENING_LOWER [eV] 0.001
  BROADENING_UPPER [eV] 0.001

  DENSITY_CUTOFF 1.0E-12

  PRINT_DENSITY_WINDOW T
  PRINT_KINETIC_ENERGY_DENSITY T
  PRINT_REGIONAL_ENERGY_DENSITY T
  PRINT_RCP T

  STRIDE 1 1 1
&END RCP
~~~

### 21.2 表面 RCP の可視化

~~~text
&RCP ON
  ENERGY_LOWER [eV] -3.0
  ENERGY_UPPER [eV] 0.0

  BROADENING_LOWER [eV] 0.001
  BROADENING_UPPER [eV] 0.001

  DENSITY_CUTOFF 1.0E-12

  PRINT_DENSITY_WINDOW F
  PRINT_KINETIC_ENERGY_DENSITY T
  PRINT_REGIONAL_ENERGY_DENSITY F
  PRINT_RCP T

  MPI_IO F
  STRIDE 1 1 1
&END RCP
~~~

### 21.3 Near-metallic / metallic 系

SCF:

~~~text
&SMEAR
  METHOD FERMI_DIRAC
  ELECTRONIC_TEMPERATURE [K] 300
&END SMEAR
~~~

RCP window broadening の例:

~~~text
BROADENING_LOWER [eV] 0.01
BROADENING_UPPER [eV] 0.01
~~~

SCF の electronic temperature と、RCP energy-window の broadening は別のパラメータである。

---

## 22. 検証済みサンプル

### benzene

~~~text
DEV_TREE/local-sham/rcp_validation/benzene/
~~~

特徴:

- closed-shell
- nonperiodic
- \([-3,0]\) eV window

### C2H5

~~~text
DEV_TREE/local-sham/rcp_validation/c2h5/
~~~

特徴:

- UKS doublet
- nonperiodic
- \([-2,0]\) eV window

### 5chain

~~~text
DEV_TREE/local-sham/rcp_validation/5chain/
~~~

特徴:

- 210 C atoms
- periodic
- Γ point
- near-metallic
- Fermi–Dirac 300 K
- \([-1,0]\) eV window
- final cube は MPI_IO F

---

## 23. Quick checklist

新しい系で RCP を計算するときは次を確認する。

1. GPW であり、Γ点または検証済みの k-point sampling を使う。
2. RKS または collinear UKS である。
3. SCF が十分収束している。
4. energy window の物理的意味を決める。
5. ENERGY_LOWER < ENERGY_UPPER である。
6. BROADENING_LOWER / UPPER は正の値にする。
7. RCP が必要なら PRINT_RCP T。
8. interface 可視化には PRINT_KINETIC_ENERGY_DENSITY T。
9. 大規模 cube が止まる場合は MPI_IO F。
10. 定量比較では STRIDE 1 1 1。
11. window electron count の3種類の値が一致することを確認する。
12. Integral(T_e) と E_kin が一致することを確認する。
13. vacuum / node の RCP は DENSITY_CUTOFF の影響を考慮する。
14. cube の単位を取り違えない。

---

## 24. 関連ソースコード

RCP input section:

~~~text
src/input_cp2k_print_dft.F
~~~

RCP driver:

~~~text
src/qs_energy_window.F
~~~

SCF 後の RCP 呼び出し:

~~~text
src/qs_scf_post_gpw.F
~~~

cube writer:

~~~text
src/pw/realspace_grid_cube.F
~~~

validation data:

~~~text
local-sham/rcp_validation/
~~~

---

## 25. 参考文献

RCP の理論的背景については次を参照する。

M. Fukuda et al., Regional chemical potential analysis for material surfaces,
J. Chem. Phys. 164, 084123 (2026).
DOI: 10.1063/5.0288934
