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

  ``use_ibrav`` (Quantum ESPRESSO)

    結晶構造の入力に Quantum ESPRESSO の ``ibrav`` パラメータを利用します。 ``true`` の場合、格子のとり方を Quantum ESPRESSO の convention に合うように変換します。入力ファイルにはあわせて格子に関するパラメータ ``a``, ``b``, ``c``, ``cosab``, ``cosac``, ``cosbc`` が(必要に応じて)書き出されます。
    (デフォルト値: ``false``)

  ``tolerance``

    ``use_ibrav = true`` の場合に、再構成した Structure データと元データとの一致を評価する際の許容度を指定します。
    (デフォルト値: 0.01)

  ``supercell``

    supercell を設定する場合に supercell のサイズを [:math:`n_x`, :math:`n_y`, :math:`n_z`] で指定します。
    (デフォルト値: なし)

  ``sqs_transformation`` ブロック

    Special Quasi-random Structure (SQS) 法を用いて合金や混晶の元素配置を求めます。この機能を利用するには `Alloy Theoretic Automated Toolkit (ATAT) <https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/>`_ が必要です。

    以下のパラメータを指定できます。パラメータの詳細は `pymatgen のマニュアル <https://pymatgen.org/pymatgen.transformations.html#pymatgen.transformations.advanced_transformations.SQSTransformation>`_ を参照してください。

      ``enable``

        SQS transformation を適用します。(デフォルト値: False)

      ``scaling``

        スーパーセルのサイズを整数値または整数のリストで指定します。整数値の場合、スーパーセル内の単位胞の数です。整数のリストの場合は、a, b, c 軸に沿った単位胞の数です。指定しない場合はエラーとして SQS transformation を適用しません。

      ``cluster_size_and_shell``

        クラスター法のパラメータをセットします。パラメータの詳細は pymatgen のマニュアルを参照してください。

      ``option``

        pymatgen の SQSTransformation のパラメータをテーブル形式で指定します。以下のパラメータが指定できます。デフォルト値は pymatgen に準拠します。
	``search_time``, ``directory``, ``instances``, ``temperature``, ``wr``, ``wn``, ``wd``, ``tol``, ``best_only``, ``remove_duplicate_structures``, ``reduction_algo``

	``search_time`` は探索の打ち切り時間を分単位で指定します。デフォルト値は 1分です。

	``instances`` は ATAT の mcsqs プログラムを並列実行する際の並列度を指定します。デフォルト値は 1 です。 ``all`` を指定すると実行環境のコア数と同数になります。 ``env`` を指定すると環境変数 ``OMP_NUM_THREADS`` の値を使用します。環境変数がセットされていない場合は 1 になります。

      ``output_cif``

        SQS transformation 後の構造データを CIF 形式で出力します。出力ファイル名を ``output_cif`` に指定します。

optional
--------
第一原理計算プログラムごとに必要な global な設定を行います。記述する内容は以下の各節に記述します。
    
tasks
-----
入力ファイルの内容を記述します。複数の入力ファイルに対応するため、 ``tasks`` には各入力ファイルごとのブロックからなるリスト形式をとります。各ブロックに記述する項目は以下のとおりです。

  ``mode`` (Quantum ESPRESSO)

    入力ファイルの計算内容を指定します。
    現時点では Quantum ESPRESSO の pw.x 向けに ``scf`` と ``nscf`` に対応しています。対応していない mode については、 ``content`` の内容をそのまま出力します。

  ``output_file`` (Quantum ESPRESSO)

    出力ファイル名を指定します。
    
  ``output_dir``

    出力先のディレクトリ名を指定します。デフォルト値はカレントディレクトリです。

  ``content``

    出力内容を指定します。Quantum ESPRESSO の場合は ``namelist`` ブロックに namelist データ (``&system`` や ``&control`` など) を記述し、 ``K_POINTS`` 等の card データを個別のブロックとして記述します。card データはパラメータをとるものがあります。

  ``template`` (Quantum ESPRESSO)

  ``template_dir`` (VASP)

    出力内容のテンプレートファイルを指定します。指定がない場合はテンプレートを利用しません。
    このファイルの内容を ``content`` に追加します。同じデータがある場合は後者を優先します。

パラメータセット指定
--------------------
入力パラメータには値のリストや範囲を指定することができ、値の組み合わせごとに個別のディレクトリを作成して入力ファイルを生成します。パラメータセットの指定には特別な構文 ``${...}`` を用います。
指定方法は次の通りです:

- リスト ``${[ A, B, ... ]}``

  パラメータセットを Python のリストとして列挙します。各項目はスカラー値やリストを指定できます。

- 範囲指定 ``${range(N)}``, ``${range(start, end, step)}``

  パラメータの範囲を指定します。それぞれ 0〜N-1, start〜end を step 刻み (step を省略した場合は 1) です。 ``N``, ``start``, ``end``, ``step`` は int または float です。

Quantum ESPRESSO 向けパラメータ
===============================

``optional`` セクションおよび ``tasks`` セクションの ``content`` に記載する内容について、Quantum ESPRESSO 固有の内容を記述します。
現時点では pw.x の ``scf`` および ``nscf`` に対応しています。

optionalセクション
------------------

  ``pp_file``

    元素種と擬ポテンシャルを対応付けるCSV形式のインデックスファイルを指定します。
    ファイルの書式は、元素種、擬ポテンシャルファイルのタイプ、nexclude、orbitals です。例:

    .. code-block::

      Fe,pbe-spn-rrkjus_psl.0.2.1,4,spd

    擬ポテンシャルファイルのファイル名は Fe.pbe-spn-rrkjus_psl.0.2.1.UPF に対応します。

  ``cutoff_file``

    擬ポテンシャルファイルとカットオフを対応付けるCSV形式のインデックスファイルを指定します。
    ファイルの書式は、擬ポテンシャルファイル、ecutwfcの値、ecutrhoの値 です。

  ``pseudo_dir``

    擬ポテンシャルファイルを格納するディレクトリ名を指定します。カットオフの値を擬ポテンシャルファイルから取得する場合に使用します。Quantum ESPRESSO の ``pseudo_dir`` パラメータとは独立に指定します。


content
--------

  **namelist**

  - ``structure`` セクションの ``use_ibrav`` パラメータに応じて、 ``&system`` の格子情報の指定が上書きされます。

    - ``use_ibrav = false``:
      ``ibrav`` は 0 にセットされます。また、格子パラメータに関する ``a``, ``b``, ``c``, ``cosab``, ``cosac``, ``cosbc``, ``celldm`` は削除されます。

    - ``use_ibrav = true``:
      ``ibrav`` は結晶構造データから取得された Bravais 格子のインデックスがセットされます。また、Structure データは基本格子のとり方など Quantum ESPRESSO の convention に合わせて再構成されます。

  - ``nat`` (原子数) および ``ntyp`` (元素種の数)は結晶構造データから取得される値で上書きされます。

  - カットオフ ``ecutwfc`` および ``ecutrho`` の情報は、パラメータの値が空欄の場合は擬ポテンシャルファイルから取得します。

  **CELL_PARAMETERS**

  - ``structure`` セクションの ``use_ibrav`` パラメータが true の場合は出力されません。false の場合は格子ベクトルが出力されます。単位は angstrom です。

  - 格子ベクトルの情報は結晶構造データから取得されます。 ``data`` フィールドに 3x3 の行列を直接指定した場合はその値が用いられます。

  **ATOMIC_SPECIES**

  - 原子種・原子量・擬ポテンシャルファイル名のリストを出力します。

  - 原子種の情報は結晶構造データから取得されます。擬ポテンシャルのファイル名は ``pp_list`` で指定するCSV形式のインデックスファイルを参照します。

  - ``data`` フィールドに必要なデータを指定した場合はその値が用いられます。

  **ATOMIC_POSITIONS**

  - 原子種と原子座標(fractional coordinate)のリストを出力します。

  - ``ignore_species`` に原子種または原子種のリストを指定した場合、その原子種については ``if_pos`` の値が 0 にセットされます。MDや構造最適化の際に使われます。

  - ``data`` フィールドに必要なデータを指定した場合はその値が用いられます。

  **K_POINTS**

  - k点の情報を出力します。 ``option`` に出力タイプを指定します。

    - ``gamma``: :math:`\Gamma` 点を用います。

    - ``crystal``: メッシュ状の k点のリストを出力します。メッシュの指定は ``grid`` パラメータまたは ``vol_density`` や ``k_resolution`` から導出される値が用いられます。

    - ``automatic``: k点のメッシュを指定します。メッシュの指定は ``grid`` パラメータまたは ``vol_density`` や ``k_resolution`` から導出される値が用いられます。 シフトの指定は ``kshifts`` パラメータを参照します。
    
  - メッシュの指定は以下の順序で決定されます。

    - ``grid`` パラメータの指定。grid の値は :math:`n_x, n_y, n_z` の配列またはスカラー値 :math:`n` です。後者の場合は :math:`n_x = n_y = n_z = n` と仮定します。
    - ``vol_density`` パラメータから自動導出。
    - ``k_resolution`` パラメータから自動導出。``k_resolution`` のデフォルトは 0.15 です。

  - ``data`` フィールドに必要なデータを指定した場合はその値が用いられます。


VASP 向けパラメータ
===============================

``optional`` セクションおよび ``tasks`` セクションの ``content`` に記載する内容について、VASP 固有の内容を記述します。

optional
--------

擬ポテンシャルのタイプや格納場所を指定します。

pymatgen では、擬ポテンシャルファイルを
``PMG_VASP_PSP_DIR``/*functional*/POTCAR. *element* (.gz) または
``PMG_VASP_PSP_DIR``/*functional*/ *element* /POTCAR から取得します。
``PMG_VASP_PSP_DIR`` はディレクトリの指定で、設定ファイル ``~/.config/.pmgrc.yaml`` に記載するか、同名の環境変数に指定します。また、 *functional* は擬ポテンシャルのタイプで、 ``POT_GGA_PAW_PBE`` や ``POT_LDA_PAW`` などが決められています。

  ``pseudo_functional``

    擬ポテンシャルのタイプを指定します。タイプの指定と上記の *functional* の対応は pymatgen 内のテーブルに定義され、 ``PBE`` → ``POT_GGA_PAW_PBE``,　 ``LDA`` → ``POT_LDA_PAW`` などのようになっています。

以下の ``pseudo_dir`` を指定した場合は pymatgen の流儀を無視して擬ポテンシャルの格納ディレクトリを探します。
    
  ``psuedo_dir``

    擬ポテンシャルの格納ディレクトリを指定します。擬ポテンシャルファイルのファイル名は ``pseudo_dir``/POTCAR. *element* (.gz) または ``pseudo_dir``/*element*/POTCAR です。


tasks
-----

テンプレートファイルは、 ``template_dir`` で指定するディレクトリ内に ``INCAR``, ``KPOINTS``, ``POSCAR``, ``POTCAR`` ファイルを配置します。ファイルがない項目は無視されます。

content
-------

  **incar**

  - INCAR ファイルに記述するパラメータを列挙します。

  **kpoints**

  - ``type``

    KPOINTS の指定方法を記述します。以下の値に対応しています。タイプによりパラメータが指定可能なものがあります。詳細は pymatgen.io.vasp のマニュアルを参照してください。

    - ``automatic``

      parameter: ``grid``

    - ``gamma_automatic``

      parameter: ``grid``, ``shift``

    - ``monkhorst_automatic``

      parameter: ``grid``, ``shift``

    - ``automatic_density``

      parameter: ``kppa``, ``force_gamma``

    - ``automatic_gamma_density``

      parameter: ``grid_density``

    - ``automatic_density_by_vol``

      parameter: ``grid_density``, ``force_gamma``

    - ``automatic_density_by_lengths``

      parameter: ``length_density``, ``force_gamma``

    - ``automatic_linemode``

      parameter: ``division``, ``path_type`` (HighSymmKpath の path_type に対応)


OpenMX 向けパラメータ
===============================

``optional`` セクションおよび ``tasks`` セクションの ``content`` に記載する内容について、OpenMX 固有の内容を記述します。

optional
--------

  ``data_path``

    擬原子軌道および擬ポテンシャルのファイルを格納するディレクトリを指定します。入力ファイルの ``DATA.PATH`` パラメータに対応します。

content
--------

  ``precision``

    擬原子軌道を OpenMXマニュアル 10.6 章の Table 1, 2 にしたがって選択します。 ``quick``, ``standard``, ``precise`` のいずれかの値を取ります。デフォルト値は ``quick`` です。

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
