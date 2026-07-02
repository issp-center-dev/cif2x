.. _sec-cif2x-tutorial:

チュートリアル
================================================================

第一原理計算入力ファイル生成ツール cif2x を使うには、入力パラメータファイルと結晶構造データおよび擬ポテンシャルファイルを用意した後、プログラム cif2x を実行します。現在は Quantum ESPRESSO, VASP, OpenMX, AkaiKKR の入力ファイル生成に対応しています。以下では、 ``docs/tutorial/cif2x`` ディレクトリにある Quantum ESPRESSO 向けサンプルを例にチュートリアルを実施します。

入力パラメータファイルを作成する
----------------------------------------------------------------

入力パラメータファイルには、第一原理計算プログラムに与える入力ファイルの内容を記述します。

以下に入力パラメータファイルのサンプルを記載します。このファイルは YAML形式のテキストファイルで、結晶構造データに対するオプションの指定や、出力する第一原理計算入力ファイルの内容を記述します。仕様の詳細については :ref:`ファイルフォーマット <sec-cif2x-fileformat>` の章を参照してください。

YAMLフォーマットでは、 ``keyword: value`` の辞書形式でパラメータを記述します。 ``value`` には数値や文字列などのスカラー値や、複数の値を ``[ ]`` または箇条書きの形式で列挙するリスト型、または辞書型を入れ子にすることができます。

.. literalinclude:: ../../../../tutorial/cif2x/input.yaml
   :language: yaml


入力パラメータファイルは ``structure``, ``optional``, ``tasks`` のブロックから構成されます。
``structure`` は結晶構造データに関するオプションを指定します。
``optional`` は擬ポテンシャルに関する global な設定などを行います。

``tasks`` は出力する第一原理計算入力ファイルの内容を指定します。一連の計算に対応して複数のファイルを生成できるよう、tasks は配列の値を取ります。
各出力について、計算内容は ``mode`` パラメータで指定します。SCF計算の ``scf`` や NSCF計算の ``nscf`` に対応するほか、一般的な出力を行う任意の出力モードを指定できます。
ファイルは ``output_dir`` および ``output_file`` で指定するファイルに書き出されます。

出力内容は ``content`` に記載します。
Quantum ESPRESSO の入力ファイルは、 ``&keyword`` で始まる Fortran90 の namelist 形式と、 ``K_POINTS`` などのキーワードで始まり空行で分割される cards と呼ばれるブロックからなります。 ``content`` には namelist と cards を入れ子の辞書形式で指定します。
いくつかの例外を除いて、指定された内容が基本的にはそのまま入力ファイルに書き出されます。値が空欄のキーワードは、結晶構造データなどから求めた値が代入されます。

また、template として入力ファイルの雛形を指定することもできます。 ``template`` に指定したファイルの内容と ``content`` のデータを合わせたものを入力データとして扱います。同じキーワードのデータが存在する場合は ``content`` の指定が優先されます。従って、既存の入力ファイルを元に必要な箇所を入力パラメータファイルで上書きする使い方が可能です。上記の例では次のファイル(``scf.in_tmpl``)を template として取り込み、カットオフと CELL_PARAMETER, ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINT の箇所を結晶構造等から決めます。 ``ecutwfc`` と ``ecutrho`` が空欄で上書きされていることに留意してください。

.. literalinclude:: ../../../../tutorial/cif2x/scf.in_tmpl
   :language: fortran

	      
第一原理計算入力ファイルを生成する
----------------------------------------------------------------

入力パラメータファイル(input.yaml)と結晶構造データ(Co3SnS2_nosym.cif)を入力として cif2x を実行します。

.. code-block:: bash

  $ cif2x -t QE input.yaml Co3SnS2_nosym.cif

実行にあたり、擬ポテンシャルに関する以下の準備が必要です。

- **擬ポテンシャルファイル (.UPF) の配置**: 計算に用いる Quantum ESPRESSO の擬ポテンシャルファイル (``.UPF`` 形式) を、入力パラメータファイルの ``optional.pseudo_dir`` で指定したディレクトリ (この例では ``./pseudo``) に配置します。擬ポテンシャルファイル自体はリポジトリには含まれていないため、PSlibrary や Quantum ESPRESSO の公式擬ポテンシャル配布ページなどから、使用する擬ポテンシャルライブラリに対応するファイルを入手してください。

- **インデックスファイル (CSV) の配置**: 元素種と擬ポテンシャル名を対応付けるインデックスファイルを ``optional.pp_file`` (この例では ``./pseudo/pp_psl_pbe_rrkjus.csv``) に指定します。このサンプルで用いる PSlibrary (PBE, RRKJUS) 用のインデックスファイルは ``docs/tutorial/cif2x/pseudo/pp_psl_pbe_rrkjus.csv`` として同梱されています。別の擬ポテンシャルライブラリを使う場合は、これにならって作成してください。

なお、 ``optional.pseudo_dir`` は cif2x が ``.UPF`` ファイルを探してカットオフ等の情報を取得するためのディレクトリで、生成される入力ファイル内の Quantum ESPRESSO の ``control.pseudo_dir`` (計算実行時に pw.x が擬ポテンシャルを参照するパス) とは独立です。後者は ``content`` の ``namelist`` で指定し、空欄にしておくと cif2x が ``optional.pseudo_dir`` の値で補完します (この例では ``./pseudo`` が書き込まれます)。

``cif2x`` を実行すると Quantum ESPRESSO用の入力ファイルが生成され出力されます。出力先は入力パラメータファイル内のパラメータで指定するディレクトリ(``output_dir``)およびファイル(``output_file``)です。この例では ``./scf/scf.in`` に SCF計算用の入力ファイルが書き出されます。生成されたファイルの内容は次のようになります。

.. literalinclude:: ../../../../tutorial/cif2x/scf/scf.in
   :language: fortran
   :force:

カットオフ ``ecutwfc``, ``ecutrho`` の値が擬ポテンシャルファイルから取得され、 ``CELL_PARAMETERS``, ``ATOMIC_POSITIONS`` および ``ATOMIC_SPECIES`` の元素種が結晶構造データから決まり、 ``ATOMIC_SPECIES`` の擬ポテンシャルファイル名が ``optional.pp_file`` の対応付けから決まり、 ``K_POINTS`` が ``content.K_POINTS`` の設定(この例では ``grid: [8,8,8]``)から生成されていることが確認できます。

VASP, OpenMX, AkaiKKR 向けのサンプルは、リポジトリの ``sample/cif2x/`` ディレクトリ以下に用意されています。

パラメータセットを指定する
----------------------------------------------------------------

入力パラメータ内の値をいくつか変えながら一連の入力ファイルを生成したいことがあります。例えばカットオフの値やk点の数を変えて収束性を評価するなどの場合です。入力パラメータには値のリストや範囲を指定することができ、値の組み合わせごとに個別のディレクトリを作成して入力ファイルを生成します。
パラメータセットの指定は特別な構文 ``${...}`` を用います。

.. code-block:: yaml

   content:
     K_POINTS:
       option: automatic
       grid:   ${ [ [4,4,4], [8,8,8], [12,12,12] ] }

例えば上記のように ``K_POINTS`` を指定すると、 ``grid`` の値が ``[4,4,4]``, ``[8,8,8]``, ``[12,12,12]`` の入力ファイルがそれぞれ ``4x4x4/``, ``8x8x8/``, ``12x12x12/`` サブディレクトリ内に作成されます。

同じ ``${...}`` 構文は Quantum ESPRESSO だけでなく全ターゲットで利用でき、出力
サブディレクトリ名はスイープ値から生成されます。以下の各 ``content:`` ブロックは、
前述の例と同様に ``tasks:`` の各エントリ配下に置かれます。スイープ対象のキーは
``content`` 内でそのパラメータが定義される階層に記述します。VASP では ``incar``
(または ``kpoints``/``poscar``/``potcar``)の下、OpenMX と AkaiKKR では ``content``
直下のフラットなキーです。

VASP(カットオフ収束):

.. code-block:: yaml

   content:
     incar:
       ENCUT: ${ [400, 600, 800] }

により ``400/``, ``600/``, ``800/`` に入力ファイルが生成されます。

OpenMX(エネルギーカットオフ収束):

.. code-block:: yaml

   content:
     scf.energycutoff: ${ [150, 200, 250] }

により ``150/``, ``200/``, ``250/`` に入力ファイルが生成されます。

AkaiKKR(ブリルアンゾーン品質の収束):

.. code-block:: yaml

   content:
     bzqlty: ${ [12, 16, 20] }

により ``12/``, ``16/``, ``20/`` に入力ファイルが生成されます。

RESPACK ワークフロー
----------------------------------------------------------------

cif2x は、cRPA 計算パッケージ RESPACK に至る一連のワークフロー (Quantum ESPRESSO による SCF/NSCF 計算 → ``qe2respack`` による変換 → RESPACK 実行) のための入力ファイルをまとめて生成できます。cif2x の役割は入力ファイルの生成までで、計算の実行そのものは行いません。以下では ``docs/tutorial/cif2x/respack`` ディレクトリにある SrVO3 (立方晶ペロブスカイト) のサンプルを例に、3 つの入力ファイル (``scf/scf.in``, ``nscf/nscf.in``, ``input.in``) を生成します。

入力パラメータファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../../tutorial/cif2x/respack/input.yaml
   :language: yaml

``tasks`` には ``scf``, ``nscf``, ``respack`` の 3 つのタスクを記述します。``scf`` / ``nscf`` タスクは Quantum ESPRESSO 向けと同じ書式 (``content.namelist`` と cards) で記述します。RESPACK 固有の注意点は次のとおりです。

- ``structure.use_primitive: true`` と ``use_ibrav: false`` は必須です。RESPACK 用の高対称 k 経路 (pymatgen の ``HighSymmKpath``) が標準化された primitive cell を前提とするためです。
- ``nscf`` タスクでは ``system`` に ``nosym: true`` / ``noinv: true`` を指定します。``qe2respack`` が対称性で畳まれていない一様 k メッシュを要求するためで、``K_POINTS`` も ``option: crystal`` により全 k 点を列挙する形式にします。
- ``nbnd`` は Wannier 関数の構成と cRPA 遮蔽計算に必要な空バンド数を確保する物理パラメータで、ユーザーが指定します (この例では 60)。

RESPACK テンプレート
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``respack`` タスクは、``template`` に指定した RESPACK の namelist 雛形と結晶構造から ``input.in`` を生成します。テンプレートは ``&param_*`` namelist のみで構成します (namelist の外側に内容を書くことはできません)。

.. literalinclude:: ../../../../tutorial/cif2x/respack/respack.in_tmpl
   :language: fortran

物理量はユーザーがテンプレートで与えます: ``N_wannier = 3`` (V の t2g 軌道)、エネルギー窓 ``Lower_energy_window`` / ``Upper_energy_window``、補間 k メッシュ ``dense``、および ``&param_chiqw`` の cRPA 設定 (``flg_cRPA = 1``) です。一方、``&param_interpolation`` の高対称 k 経路 (``N_sym_points`` と座標ブロック) は cif2x が結晶構造から自動生成し、Wannier 関数の初期推定は SCDM (``N_initial_guess = 0``) に固定されます。

擬ポテンシャルの準備
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RESPACK は Kohn-Sham 波動関数を直接扱うため、ノルム保存 (ONCV) 擬ポテンシャルを使用します。この例では SG15 ONCV (PBE, v1.0) を使用します。cif2x は擬ポテンシャルファイル名を ``<元素>.<名前>.UPF`` の形式で組み立てるため、ダウンロード後にリネームして ``optional.pseudo_dir`` (この例では ``./pseudo``) に配置します。

.. code-block:: bash

  $ mkdir -p pseudo
  $ cd pseudo
  $ for el in Sr V O; do
  >   curl -LO "http://www.quantum-simulation.org/potentials/sg15_oncv/upf/${el}_ONCV_PBE-1.0.upf"
  >   mv ${el}_ONCV_PBE-1.0.upf ${el}.ONCV_PBE-1.0.UPF
  > done

元素種と擬ポテンシャル名の対応付けは ``optional.pp_file`` (``pp.csv``) で行います。

.. literalinclude:: ../../../../tutorial/cif2x/respack/pp.csv

ONCV の UPF ヘッダーにはカットオフの推奨値が含まれないため、この例ではカットオフを ``optional.cutoff_file`` (``cutoff.csv``) で与えます (``ecutwfc`` = 80 Ry, ``ecutrho`` = 320 Ry。値は目安であり、実際の計算では収束を確認してください)。

.. literalinclude:: ../../../../tutorial/cif2x/respack/cutoff.csv

入力ファイルを生成する
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

  $ cif2x -t respack input.yaml SrVO3.cif

SCF 計算用の入力ファイルが ``scf/scf.in`` に書き出されます。

.. literalinclude:: ../../../../tutorial/cif2x/respack/scf/scf.in
   :language: fortran
   :force:

NSCF 計算用の入力ファイル ``nscf/nscf.in`` では、``calculation = 'nscf'`` に加えて ``nosym``, ``noinv``, ``nbnd`` が設定され、``K_POINTS crystal`` に全 k 点が列挙されます (紙面の都合で k 点リストは省略します)。

RESPACK 用の制御ファイルは ``input.in`` に書き出されます。``&param_interpolation`` の namelist 終端 (``/``) の直後に、高対称 k 経路の座標ブロックが自動生成されている点に注意してください。

.. literalinclude:: ../../../../tutorial/cif2x/respack/input.in
   :language: fortran

次のステップ: RESPACK の実行
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

生成した入力ファイルを使った計算は次の順で実行します (コマンドのみ示します。実行方法・並列化の詳細は Quantum ESPRESSO および RESPACK のドキュメントを参照してください)。

.. code-block:: bash

  $ pw.x -in scf/scf.in > scf.out                # SCF 計算
  $ pw.x -in nscf/nscf.in > nscf.out             # NSCF 計算
  $ qe2respack.py work/pwscf.save                # QE 出力を RESPACK 形式に変換
  $ calc_wannier < input.in > log.wannier        # Wannier 関数の構成
  $ calc_chiqw < input.in > log.chiqw            # cRPA 分極関数
  $ calc_w3d < input.in > log.w3d                # 有効クーロン相互作用 W
  $ calc_j3d < input.in > log.j3d                # 有効交換相互作用 J
