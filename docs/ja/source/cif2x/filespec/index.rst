.. _sec-cif2x-fileformat:

======================
 ファイルフォーマット
======================

入力パラメータファイル
======================

入力パラメータファイルでは、cif2x で第一原理計算入力ファイルを生成するための設定情報を YAML形式で記述します。本ファイルは以下の部分から構成されます。

  1. structureセクション: 結晶構造データの扱いについてのオプションを記述します。

  2. optionalセクション: 擬ポテンシャルファイルの指定や、YAMLの参照機能を利用する場合のシンボル定義を行います。

  3. tasksセクション: 入力ファイルの内容を記述します。

     
structure
---------

  ``use_ibrav`` (デフォルト値: ``false``)

    結晶構造の入力に Quantum ESPRESSO の ``ibrav`` パラメータを利用します。 ``true`` の場合、格子のとり方を Quantum ESPRESSO の convention に合うように変換します。入力ファイルにはあわせて格子に関するパラメータ ``a``, ``b``, ``c``, ``cosab``, ``cosac``, ``cosbc`` が(必要に応じて)書き出されます。再構成した構造と元の構造との一致は固定の許容度 ``0.01`` で評価され、下記の ``tolerance`` パラメータはこの評価には影響しません。

  ``use_primitive`` (デフォルト値: ``false``)

    ``true`` の場合、入力ファイルを生成する前に結晶構造を基本格子 (primitive cell) に変換します (pymatgen の ``get_primitive_structure`` を使用)。高対称 k 経路が標準基本格子で定義されるため、 ``respack`` ターゲットおよび Quantum ESPRESSO の ``bands`` モードでは必須です。

  ``tolerance`` (デフォルト値: 0.01)

    ``use_primitive = true`` の場合に、基本格子への変換で用いる許容度を指定します。それ以外の場合は影響しません。

  ``supercell`` (デフォルト値: なし)

    supercell を設定する場合に supercell のサイズを [:math:`n_x`, :math:`n_y`, :math:`n_z`] で指定します。例えば各軸方向に 2倍の supercell を作成する場合は次のように指定します。

    .. code-block:: yaml

      structure:
        supercell: [2, 2, 2]

  ``cell_shape`` (分子構造データの場合は必須)

    結晶構造ではなく孤立分子の構造データ (XYZ ファイルなど) を入力する場合に、分子を配置する箱の大きさを指定します。分子は箱の中心に置かれます。スカラー値 ``L`` (一辺 ``L`` angstrom の立方体)、リスト ``[Lx, Ly, Lz]`` (各辺の長さ、単位は angstrom)、または ``auto`` を指定します。 ``auto`` の場合は仮に一辺 10 angstrom の立方体が使われ、OpenMX では ``Atoms.UnitVectors`` ブロックを出力せず、セルの決定を OpenMX に任せます。分子構造データに対して ``cell_shape`` を指定しない場合はエラーで停止します。

  ``tol_deg`` (デフォルト値: 5.0)

    結晶構造データが磁気モーメントを持つ場合 (mCIF ファイルなど) に、モーメントが共線 (collinear) かどうかを判定する角度の許容度を度単位で指定します。共線でない場合、Quantum ESPRESSO では ``nspin = 2`` が指定されていても非共線計算に切り替えられます。

optional
--------
第一原理計算プログラムごとに必要な global な設定を行います。記述する内容は以下の各節に記述します。各タスクにも同じキーを持つ ``optional`` ブロックを置くことができ、タスク側の指定が global の指定より優先されます。例えばあるタスクだけ別の ``pseudo_dir`` を使うことができます。値がマッピングの項目 (``vps_variant`` など) は全体を置き換えるのではなくキー単位で統合されます。

tasks
-----
入力ファイルの内容を記述します。複数の入力ファイルに対応するため、 ``tasks`` には各入力ファイルごとのブロックからなるリスト形式をとります。各ブロックに記述する項目は以下のとおりです。

  ``mode`` (Quantum ESPRESSO, RESPACK)

    入力ファイルの計算内容を指定します。
    Quantum ESPRESSO の pw.x 向けには ``scf``, ``nscf``, ``relax``, ``vc-relax``, ``bands`` に対応しており、値は ``&control`` の ``calculation`` パラメータに書き出されます。対応していない mode については、 ``content`` の内容をそのまま出力します。 ``respack`` ターゲットでは、 ``scf``, ``nscf``, ``bands`` のタスクは Quantum ESPRESSO の入力として生成され、 ``respack`` のタスクが RESPACK の入力を生成します (:ref:`sec-cif2x-filespec-respack` 参照)。

  ``output_file`` (Quantum ESPRESSO, OpenMX, AkaiKKR, RESPACK では必須)

    出力ファイル名を指定します。VASP ではファイル名が固定 (``INCAR``, ``KPOINTS``, ``POSCAR``, ``POTCAR``) のため使用しません。

  ``output_dir``

    出力先のディレクトリ名を指定します。デフォルト値はカレントディレクトリです。

  ``content``

    出力内容を指定します。Quantum ESPRESSO の場合は ``namelist`` ブロックに namelist データ (``&system`` や ``&control`` など) を記述し、 ``K_POINTS`` 等の card データを個別のブロックとして記述します。card データはパラメータをとるものがあります。

  ``template`` (Quantum ESPRESSO, OpenMX, AkaiKKR, RESPACK)

  ``template_dir`` (VASP)

    出力内容のテンプレートファイルを指定します。指定がない場合はテンプレートを利用しません。
    このファイルの内容を ``content`` に追加します。同じデータがある場合は後者を優先します。

  ``precision`` (OpenMX)

    擬原子軌道を OpenMXマニュアル 10.6 章の Table 1, 2 にしたがって選択します。 ``quick``, ``standard``, ``precise`` のいずれかの値を取ります。デフォルト値は ``quick`` です。 ``content`` の中ではなく、タスク自体のパラメータとして記述します。

  ``optional``

    global な ``optional`` セクション (上記) をタスク単位で上書きします。キーは同じで、ここに記述した項目はこのタスクに限って global の指定より優先されます。

パラメータセット指定
--------------------
入力パラメータには値のリストや範囲を指定することができ、値の組み合わせごとに個別のディレクトリを作成して入力ファイルを生成します。パラメータセットの指定には特別な構文 ``${...}`` を用います。
指定方法は次の通りです:

- リスト ``${[ A, B, ... ]}``

  パラメータセットを Python のリストとして列挙します。各項目はスカラー値やリストを指定できます。引用符のない語は文字列として扱われます (``${[ pbe, lda ]}``)。

- 範囲指定 ``${range(N)}``, ``${range(start, end, step)}``

  パラメータの範囲を指定します。前者は 0〜N-1、後者は start から step 刻みで end の **直前まで** (end は含みません。Python の ``range`` と同じ) です。step を省略した場合は 1 です。 ``N``, ``start``, ``end``, ``step`` は int または float で、float の刻みでは値が (小数第 10 位で) 丸められるため、 ``${range(0.1, 0.4, 0.1)}`` はちょうど ``0.1``, ``0.2``, ``0.3`` になります。

以下の規則があります。

- ``${...}`` はパラメータの値全体でなければなりません (``ecutwfc: ${range(30, 60, 10)}``)。文字列の一部に ``${...}`` を含めること (``prefix: run_${[1, 2]}``) はエラーになります。

- タスク内の各 ``${...}`` は独立に展開され、すべての組み合わせについて入力ファイルが生成されます。各組み合わせのファイルは ``output_dir`` 直下に、値を出現順に ``_`` で連結した名前のサブディレクトリに書き出されます。リスト値の要素は ``x`` で連結されます (例: ``ecutwfc: 40`` と ``grid: [8, 8, 8]`` の組み合わせは ``40_8x8x8``)。

- 展開は Quantum ESPRESSO、VASP、OpenMX、AkaiKKR のタスクの ``content`` に適用されます。 ``mode: respack`` のタスクには適用されません。

Quantum ESPRESSO 向けパラメータ
===============================

``optional`` セクションおよび ``tasks`` セクションの ``content`` に記載する内容について、Quantum ESPRESSO 固有の内容を記述します。
現時点では pw.x の ``scf``, ``nscf``, ``relax``, ``vc-relax``, ``bands`` の各モードに対応しています。

optionalセクション
------------------

  ``pp_file``

    元素種と擬ポテンシャルを対応付けるCSV形式のインデックスファイルを指定します。 ``ATOMIC_SPECIES`` の生成に必須です。
    ファイルの 1行目はヘッダで、列は ``element,pseudopotential[,nexclude,orbitals][,zval]`` です。各列の意味は以下のとおりです。

    - ``element``: 元素記号。
    - ``pseudopotential``: 擬ポテンシャルファイル名のうち、元素記号と拡張子 ``.UPF`` を除いた部分(name stem)。擬ポテンシャルファイル名は ``{element}.{pseudopotential}.UPF`` として組み立てられます (スピン軌道相互作用を含む非共線計算では、生成される ``ATOMIC_SPECIES`` のファイル名が ``{element}.rel-{pseudopotential}.UPF`` のように ``rel-`` 付きになります)。
    - ``nexclude``: Wannier 関数によるバンド数の見積り (``nbnd``) で除外する内殻バンド数。namelist の ``nbnd`` を空欄にした場合のみ必要です。
    - ``orbitals``: Wannier 化の対象とする軌道を ``s``, ``p``, ``d``, ``f`` の文字で並べたもの。 ``nexclude`` とともに空欄の ``nbnd`` の計算に用いられます。
    - ``zval`` (または ``z_valence``): 擬ポテンシャルの価電子数。磁性構造データの磁気モーメントを ``starting_magnetization`` に換算する場合のみ必要です。この列がない場合は ``pseudo_dir`` 内の UPF ファイルから読み取ります (UPF v2 では ``z_valence`` 属性、UPF v1 ではヘッダの "Z valence" 行)。

    例:

    .. code-block::

      element,pseudopotential,nexclude,orbitals,zval
      Fe,pbe-spn-rrkjus_psl.0.2.1,4,spd,16

    この行は擬ポテンシャルファイル ``Fe.pbe-spn-rrkjus_psl.0.2.1.UPF`` に対応します。

  ``cutoff_file``

    擬ポテンシャルファイルとカットオフを対応付けるCSV形式のインデックスファイルを指定します。
    列は ``pseudofile,ecutwfc,ecutrho`` で、それぞれ擬ポテンシャルファイル名、 ``ecutwfc`` の値、 ``ecutrho`` の値です。``pseudofile`` 列はディレクトリを含まないファイル名 (``Fe.pbe-spn-rrkjus_psl.0.2.1.UPF``) で、 ``ATOMIC_SPECIES`` に書き出される名前と同じです。このファイルは cif2x と一緒にインストールされるコマンド ``cif2x-pp-cutoff`` で UPF ファイルから作成でき、名前はこの形式で書き出されます (``cif2x-pp-cutoff -o cutoff.csv pseudo/*.UPF`` 。詳細は ``cif2x-pp-cutoff -h``)。

  ``pseudo_dir``

    擬ポテンシャルファイルを格納するディレクトリ名を指定します。カットオフの値 (磁性構造の場合は価電子数も) を UPF ファイルから取得する場合に使用します。
    また、テンプレートや ``content`` で ``&control`` の ``pseudo_dir`` パラメータが空欄の場合は、この値が転記されます。値が指定されている場合はそのまま残ります。


content
--------

  **namelist**

  - ``structure`` セクションの ``use_ibrav`` パラメータに応じて、 ``&system`` の格子情報の指定が上書きされます。

    - ``use_ibrav = false``:
      ``ibrav`` は 0 にセットされます。また、格子パラメータに関する ``a``, ``b``, ``c``, ``cosab``, ``cosac``, ``cosbc``, ``celldm`` は削除されます。

    - ``use_ibrav = true``:
      ``ibrav`` は結晶構造データから取得された Bravais 格子のインデックスがセットされます。また、Structure データは基本格子のとり方など Quantum ESPRESSO の convention に合わせて再構成されます。

  - ``nat`` (原子数) および ``ntyp`` (元素種の数)は結晶構造データから取得される値で上書きされます。磁性構造データでは、同じ元素でもモーメントが異なる原子は別の種として数えられ、すべての種のラベルに通し番号が付きます (``Fe0``, ``Fe1``, ``O0`` のように、モーメントが 1 種類だけの元素にも付きます)。

  - pw.x のモード (``scf``, ``nscf``, ``relax``, ``vc-relax``, ``bands``) では、 ``&control`` の ``calculation`` がタスクの ``mode`` にセットされ、 ``&control`` の ``pseudo_dir`` が空欄なら ``optional.pseudo_dir`` の値が転記されます (いずれもテンプレートまたは ``content`` に ``&control`` namelist がある場合のみ)。

  - カットオフ ``ecutwfc`` および ``ecutrho`` は、値が空欄の場合 (``ecutwfc`` はキー自体がない場合、 ``ecutrho`` は両方のキーがない場合も) に解決されます。元素ごとに、まず ``cutoff_file``、次に ``pseudo_dir`` 内の UPF ファイルの ``wfc_cutoff`` / ``rho_cutoff`` 属性または "Suggested minimum cutoff" 行を参照し、元素全体での最大値を書き出します。どちらの参照にも ``pp_file`` が必要です。 ``ecutwfc`` は必ず解決できなければならず、できない場合は ``cutoff_file`` または明示的な値を求めるエラーで停止します。 ``ecutrho`` が解決できない場合はキーを出力せず、 ``pw.x`` のデフォルト (``4 * ecutwfc``) に任せます。参照して得た ``ecutrho`` がこのデフォルトより小さい場合はデフォルト値まで引き上げられますが、明示的に指定した ``ecutrho`` はそのまま書き出されます。

  - ``nbnd`` が空欄の場合は、 ``pp_file`` の ``nexclude`` と ``orbitals`` 列から ``nexclude + 2 W + ceil(W/2)`` (偶数に切り上げ) として見積もられます。 ``W`` は Wannier 化する軌道の数で、内殻バンドの上におよそ 2.5 W のバンドを取ることになります。非共線計算では ``nexclude`` と ``W`` が 2 倍になります。これらの列がない場合は ``nbnd`` を明示的に指定してください。

  - ``nspin``, ``noncolin``, ``lspinorb`` は整合するように書き直されます。共線スピン分極計算では ``nspin = 2``、非共線計算では ``noncolin = .true.`` (指定があれば ``lspinorb`` も)、非磁性計算では何も出力されません。結晶構造データが磁気モーメントを持ち、これらのキーがいずれも指定されていない場合は、モーメント自体から ``nspin = 2`` または ``noncolin`` を選びます (``tol_deg`` 参照)。モーメントが非共線なのに ``nspin = 2`` が指定されている場合は ``noncolin`` に切り替えられます。

  - 磁性構造データでは、有限のモーメントを持つ各原子種 ``i`` について ``starting_magnetization(i)`` が書き出されます。値はボーア磁子単位のモーメントを擬ポテンシャルの価電子数 (``pp_file`` の ``zval`` 列、または UPF ファイル。どちらもなければいずれかを求めるエラーになります) で割ったものです。 ``±1`` を超える値は警告付きで ``±1`` に丸められます。非共線計算では ``angle1(i)`` (z 軸からの極角) と ``angle2(i)`` (x 軸からの方位角) も度単位で書き出されます。 ``use_ibrav`` の場合、モーメントはセルとともに回転されます。テンプレートや ``content`` で特定の種に指定した値 (``starting_magnetization(1)``) は保持され、他の種は生成されます。添字なしの ``starting_magnetization`` リストを指定した場合は生成全体が行われません。 ``nspin: 1`` を明示した場合は何も生成されません。非共線モーメントに対して ``noncolin: false`` を明示した場合はその指定が優先され、モーメントは最初の有限モーメントの軸に射影されます (警告が出力されます)。

  **CELL_PARAMETERS**

  - ``structure`` セクションの ``use_ibrav`` パラメータが true の場合は出力されません。false の場合は格子ベクトルが出力されます。単位は angstrom です。

  - 格子ベクトルの情報は常に結晶構造データから取得されます。入力パラメータファイルから上書きすることはできません。

  **ATOMIC_SPECIES**

  - 原子種・原子量・擬ポテンシャルファイル名のリストを出力します。

  - 原子種の情報は結晶構造データから取得されます。擬ポテンシャルのファイル名は ``pp_file`` で指定するCSV形式のインデックスファイルを参照します。このブロックの生成には ``pp_file`` の指定が必須です。

  **ATOMIC_POSITIONS**

  - 原子種と原子座標(fractional coordinate)のリストを出力します。座標の種類は ``option`` パラメータで指定します (デフォルト値: ``crystal``)。

  - ``ignore_species`` に原子種または原子種のリストを指定した場合、その原子種については ``if_pos`` の値が 0 にセットされます。MDや構造最適化の際に使われます。

  **K_POINTS**

  - k点の情報を出力します。 ``option`` に出力タイプを指定します。

    - ``gamma``: :math:`\Gamma` 点を用います。

    - ``crystal``: メッシュ状の k点のリストを出力します。メッシュの指定は ``grid`` パラメータまたは ``vol_density`` や ``k_resolution`` から導出される値が用いられます。

    - ``automatic``: k点のメッシュを指定します。メッシュの指定は ``grid`` パラメータまたは ``vol_density`` や ``k_resolution`` から導出される値が用いられます。 シフトの指定は ``kshifts`` パラメータを参照します。

    - ``crystal_b``: 結晶の高対称 k 経路 (pymatgen の ``HighSymmKpath``) からバンド計算用の経路を生成します。高対称点間の点数は ``line_npoints`` (デフォルト値: 20)、独自の経路はラベル列のリスト ``path`` で指定します。 ``bands`` モードではこのオプションが必須です。詳細は :ref:`sec-cif2x-appendix` を参照してください。

  - メッシュの指定は以下の順序で決定されます。

    - ``grid`` パラメータの指定。grid の値は :math:`n_x, n_y, n_z` の配列またはスカラー値 :math:`n` です。後者の場合は :math:`n_x = n_y = n_z = n` と仮定します。
    - ``vol_density`` パラメータから自動導出 (逆格子体積あたりの k 点数。pymatgen の ``automatic_density_by_vol`` と同じ。分割数のみを使い、シフトは ``kshifts`` で指定します)。
    - ``k_resolution`` パラメータから自動導出。 ``k_resolution`` は k 点間隔の上限 (単位は 1/angstrom、デフォルト値 0.15) で、逆格子ベクトル :math:`\mathbf{b}_i` (因子 :math:`2\pi` を含む) から :math:`n_i = \lceil |\mathbf{b}_i| / k_\mathrm{resolution} \rceil` と決まります。

  - ``data`` フィールドが参照されるのは、バンド経路のオプション ``crystal_b`` と ``tpiba_b`` (および未知の ``option`` を指定した場合で、この場合は記述どおりに出力されます) に限られます。これらの場合は ``data`` に記述した行がそのまま出力され、生成は行われません。それ以外のオプションおよび上記の他の card では、データは常に結晶構造から生成されます。


VASP 向けパラメータ
===============================

``optional`` セクションおよび ``tasks`` セクションの ``content`` に記載する内容について、VASP 固有の内容を記述します。

optional
--------

``optional`` セクションでは POTCAR ファイルの取得元を指定します。方法は 2 つあり、 ``pp_file`` が指定されている場合は 1 の方法が使われます。

1. POTCAR ファイルのインデックスファイル (``pp_file`` と ``pseudo_dir``)

  ``pp_file``

    ヘッダ行と ``element,pseudopotential`` の 2 列からなる CSV ファイルで、元素記号とその POTCAR ファイル (素のファイルまたは ``.gz``) のパスを対応付けます。相対パスは ``pseudo_dir`` が指定されていればそこからの相対パスです。構造に含まれるすべての元素の行が必要です。例:

    .. code-block::

      element,pseudopotential
      Na,Na_pv/POTCAR
      Cl,Cl/POTCAR

  ``pseudo_dir``

    ``pp_file`` 内のパスの基準となるディレクトリです。

2. pymatgen の流儀 (``pseudo_functional`` と ``pseudo_map``)

  ``pp_file`` がない場合は pymatgen と同じ方法で POTCAR ファイルを探します。すなわち ``PMG_VASP_PSP_DIR``/*functional*/POTCAR. *symbol* (.gz) または ``PMG_VASP_PSP_DIR``/*functional*/ *symbol* /POTCAR で、 ``PMG_VASP_PSP_DIR`` は設定ファイル ``~/.config/.pmgrc.yaml`` に記載するか同名の環境変数で指定し、 *functional* は ``POT_GGA_PAW_PBE`` や ``POT_LDA_PAW`` などのディレクトリ名です。この方法では ``PMG_VASP_PSP_DIR`` の設定が必須です。

  ``pseudo_functional``

    汎関数を指定します。 *functional* ディレクトリとの対応は pymatgen 内のテーブルで決まります (``PBE`` → ``POT_GGA_PAW_PBE``, ``LDA`` → ``POT_LDA_PAW`` など)。省略時は pymatgen のデフォルトが使われます。

  ``pseudo_map``

    ``element,symbol`` の 2 列からなる CSV ファイルで、元素に使用する POTCAR のシンボル (例: ``Na,Na_pv``) を対応付けます。記載のない元素は元素記号そのままのシンボルを使います。

POTCAR を組み立てられない場合 (``pp_file`` 自体が読めない、 ``pp_file`` に元素の行がない、POTCAR ファイルが読めない、 ``PMG_VASP_PSP_DIR`` が未設定、pymatgen が汎関数とシンボルに対応する POTCAR を見つけられない) は、POTCAR を欠いた入力一式が意図せず書き出されないよう、cif2x はエラーで停止します。計算ノード側で POTCAR を付加する場合など、意図的に POTCAR を省くには ``content`` の ``potcar`` に ``skip: true`` を指定します (下記参照)。

tasks
-----

テンプレートファイルは、 ``template_dir`` で指定するディレクトリ内に ``INCAR``, ``KPOINTS`` の名前で配置します。ファイルがない項目は無視されます。 ``POSCAR`` は常に結晶構造から、 ``POTCAR`` は ``optional`` の設定から生成されるため、 ``template_dir`` 内のこれらの名前のファイルは使用されません。

content
-------

  **incar**

  - INCAR ファイルに記述するパラメータを列挙します。キーの大文字小文字は区別しません。

  - 磁性構造データ (mCIF ファイルなど) の場合は、モーメントから ``MAGMOM`` が生成されます。共線モーメントでは ``POSCAR`` の原子順に 1 原子あたり 1 つの値 (あわせて ``ISPIN = 2``)、非共線モーメントでは 1 原子あたり 3 つの値 (あわせて ``LNONCOLLINEAR = .TRUE.``) です。ここで ``MAGMOM`` を指定した場合はその値がそのまま使われ、 ``ISPIN`` / ``LNONCOLLINEAR`` の指定も変更されません。

  **kpoints**

  - ``type``

    KPOINTS の指定方法を記述します。デフォルト値は ``gamma_automatic`` です。以下の値に対応しています。タイプによりパラメータが指定可能なものがあり、以下のパラメータ名は ``kpoints`` の下に記述するキー名です (大文字小文字は区別しません)。各方式の意味は pymatgen.io.vasp のマニュアルを参照してください。

    - ``automatic``

      parameter: ``grid`` (整数 1 つ、分割の長さ。デフォルト値: 1)

    - ``gamma_automatic``

      parameter: ``kpoints`` (リスト :math:`[n_x, n_y, n_z]` 。デフォルト値: ``[1, 1, 1]``), ``shift`` (デフォルト値: ``[0, 0, 0]``)

    - ``monkhorst_automatic``

      parameter: ``kpoints`` (リスト :math:`[n_x, n_y, n_z]` 。デフォルト値: ``[2, 2, 2]``), ``shift`` (デフォルト値: ``[0, 0, 0]``)

    - ``automatic_density``

      parameter: ``grid_density`` (逆格子原子あたりの k 点数。デフォルト値: 0.1), ``force_gamma`` (デフォルト値: ``false``)

    - ``automatic_gamma_density``

      parameter: ``grid_density`` (デフォルト値: 0.1)

    - ``automatic_density_by_vol``

      parameter: ``grid_density`` (逆格子体積あたりの k 点数。デフォルト値: 1), ``force_gamma`` (デフォルト値: ``false``)

    - ``automatic_density_by_lengths``

      parameter: ``length_density`` (3 つの長さのリスト。デフォルト値: ``[10.0, 10.0, 10.0]``), ``force_gamma`` (デフォルト値: ``false``)

    - ``automatic_linemode``

      parameter: ``division`` (デフォルト値: 1), ``path_type`` (HighSymmKpath の path_type に対応)

    例えば :math:`6 \times 6 \times 6` の :math:`\Gamma` 中心メッシュは次のように記述します。

    .. code-block:: yaml

      kpoints:
        type: gamma_automatic
        kpoints: [6, 6, 6]

  **poscar**

  - ``POSCAR`` は結晶構造から生成されます。 ``fix_species`` に元素記号またはそのリストを指定すると、該当する元素の原子の selective dynamics フラグが ``F F F``、それ以外が ``T T T`` になります。空のリストを指定すると、構造データが持つ ``selective_dynamics`` 属性が取り除かれます。

  **potcar**

  - ``skip: true`` を指定すると生成する一式から ``POTCAR`` を省きます (上記 ``optional`` セクション参照)。他のパラメータはありません。


OpenMX 向けパラメータ
===============================

``optional`` セクションおよび ``tasks`` セクションの ``content`` に記載する内容について、OpenMX 固有の内容を記述します。

optional
--------

  ``data_path``

    擬原子軌道および擬ポテンシャルのファイルを格納するディレクトリを指定します。入力ファイルの ``DATA.PATH`` パラメータに対応します。

  ``vps_variant``

    一部の元素 (Fe, Co, Ni, Cu, Zn) には OpenMX のデータベースに 2 種類の擬ポテンシャル、hard 版 (``Fe_PBE19H``、価電子数が多い) と soft 版 (``Fe_PBE19S``) があります。このパラメータで元素ごとに使用する版を、接尾辞 ``H`` / ``S`` または完全な名前で指定します (大文字小文字は区別しません)。指定がない場合は hard 版が使われます。

    .. code-block:: yaml

      optional:
        data_path: ../openmx3.9/DFT_DATA19
        vps_variant:
          Fe: S

tasks
-----

  ``precision``

    擬原子軌道を OpenMXマニュアル 10.6 章の Table 1, 2 にしたがって選択します。 ``quick``, ``standard``, ``precise`` のいずれかの値を取ります。デフォルト値は ``quick`` です。 ``content`` の中ではなく、 ``output_file`` や ``template`` と同じ階層に記述します。

    .. code-block:: yaml

      tasks:
        - output_file: openmx.dat
          template: template.in
          precision: precise
          content:
            system.name: NaCl

content
--------

``content`` には OpenMX のキーワードを記述し、そのまま出力ファイルに書き出されます (リスト値は ``<Key ... Key>`` のブロック形式で出力されます)。以下のキーワードは結晶構造データから生成され、テンプレートや ``content`` の値を上書きします: ``Species.Number``, ``Definition.of.Atomic.Species``, ``Atoms.Number``, ``Atoms.SpeciesAndCoordinates.Unit``, ``Atoms.SpeciesAndCoordinates``, ``Atoms.UnitVectors.Unit``, ``Atoms.UnitVectors`` 。 ``System.Name`` が指定されていない場合は構造ファイルのベース名が使われます。

``Atoms.SpeciesAndCoordinates`` の up / down 電子数は、非磁性構造では擬ポテンシャルの価電子数の半分ずつです。共線モーメントを持つ磁性構造データでは、原子のモーメント :math:`m` (ボーア磁子単位) から :math:`(N \pm m)/2` となり、テンプレートや ``content`` で指定がなければ ``scf.SpinPolarization`` が ``on`` にセットされます。非共線モーメントは OpenMX では未対応で、cif2x はエラーで停止します。この場合はブロックを手で記述してください。

AkaiKKR 向けパラメータ
===============================

``optional`` セクションおよび ``tasks`` セクションの ``content`` に記載する内容について、AkaiKKR 固有の内容を記述します。

optional
--------

  ``workdir``

    一時ファイルの出力先を指定します。指定しない場合は ``/tmp`` または ``TMPDIR`` 環境変数の値を利用します。


content
--------

``content`` には AkaiKKR の入力パラメータの内容を記述します。指定のない項目は空欄が出力され、AkaiKKR内部のデフォルト値が使われます。
以下のパラメータは結晶構造データから決まる値で置き換えられます。

- ``brvtyp``:
  ただし、 ``brvtyp`` に ``aux`` (を含む)値が指定されている場合は上書きされません。

- 格子パラメータ: ``a``, ``c/a``, ``b/a``, ``alpha``, ``beta``, ``gamma``, ``r1``, ``r2``, ``r3``

- タイプ情報: ``ntyp``, ``type``, ``ncmp``, ``rmt``, ``field``, ``mxl``, ``anclr``, ``conc``

- 元素種の情報: ``natm``, ``atmicx``, ``atmtyp``

なお、 ``rmt`` と ``field`` の値は、入力パラメータファイル内で ``ntyp`` 個の要素からなるリストとして指定されている場合のみ、入力パラメータファイルの値が使われます。


.. _sec-cif2x-filespec-respack:

RESPACK 向けパラメータ
===============================

``respack`` ターゲットは RESPACK ワークフローの入力一式、すなわち Quantum ESPRESSO の ``scf`` と ``nscf`` (および必要に応じて ``bands``) の入力と、それに続く RESPACK の制御ファイルを生成します。RESPACK の入力に書き出される高対称 k 経路は標準基本格子で定義されるため、 ``structure`` セクションに ``use_primitive: true`` と ``use_ibrav: false`` を指定する必要があります。

tasks
-----

各タスクは ``mode`` で区別します。

- ``scf``, ``nscf``, ``bands``: Quantum ESPRESSO ターゲットと全く同じ方法で生成されます。パラメータは上記の :ref:`sec-cif2x-fileformat` の記述に従います (``optional`` の ``pp_file``, ``cutoff_file``, ``pseudo_dir`` と、 ``content`` の ``namelist`` および card)。 ``nscf`` タスクでは ``qe2respack`` が要求する ``nosym: true`` と ``noinv: true`` を ``content.namelist.system`` に指定してください。指定がない場合 cif2x は警告を出力します。

- ``respack``: RESPACK の制御ファイルを生成します。 ``output_file`` と ``template`` が必須です。このタスクではパラメータセット指定 ``${...}`` は展開されません。

template と content (``mode: respack``)
----------------------------------------

テンプレートは namelist のみからなるファイルで、RESPACK の namelist ``&param_chiqw``, ``&param_wannier``, ``&param_interpolation``, ``&param_visualization``, ``&param_calc_int`` だけを含めることができます。これ以外の namelist や、namelist の外側のテキスト (k 点ブロックなど) はエラーになります。k 経路のブロックは cif2x が生成するためです。

``content`` には同じ namelist をマッピングとして記述します。値はテンプレートの値の上に重ねられ、同じキーについては ``content`` の値が優先されます。namelist 名とパラメータ名の大文字小文字は区別しません。

.. code-block:: yaml

  - mode: respack
    output_file: input.in
    template: respack.in_tmpl
    content:
      param_wannier:
        N_wannier: 3
      param_interpolation:
        dense: [8, 8, 8]

以下の条件が検査され、満たさない場合はエラーで停止します。

- ``&param_wannier`` の ``N_wannier`` は正の整数でなければなりません。Wannier 関数の数はユーザが決める物理量で、自動的には導出されません。

- ``&param_wannier`` の ``Lower_energy_window`` と ``Upper_energy_window`` は数値で、 ``Lower_energy_window < Upper_energy_window`` でなければなりません。

- ``&param_wannier`` の ``N_initial_guess`` は ``0`` か省略でなければなりません。初期推定は常に SCDM で、 ``N_initial_guess = 0`` が書き出されます。

- ``&param_interpolation`` の ``dense`` は整数 :math:`n` (:math:`n \times n \times n` として扱われます) または 3 つの整数のリストでなければなりません。

- ``&param_interpolation`` の ``reading_sk_format`` は ``0`` (デフォルト値) でなければなりません。

以下は cif2x が生成して出力に書き込みます: ``&param_interpolation`` の ``N_sym_points`` (経路の各区間に含まれる高対称点の数) と、 ``&param_interpolation`` の直後に続く高対称 k 点座標のブロック (pymatgen の ``HighSymmKpath`` から取得)。namelist は ``&param_chiqw``, ``&param_wannier``, ``&param_interpolation``, ``&param_visualization``, ``&param_calc_int`` の順に出力され、テンプレートまたは ``content`` に現れたものだけが出力されます。

完全な例はリポジトリの ``docs/tutorial/cif2x/respack/`` および ``sample/cif2x/respack/`` にあります。
