.. _sec-cif2x-appendix:

================================================================
拡張ガイド
================================================================

Quantum ESPRESSO の mode を追加する
----------------------------------------------------------------

Quantum ESPRESSO の計算モードへの対応を追加するには、 ``src/cif2x/qe/calc_mode.py`` の ``create_modeproc()`` 関数に mode と変換クラスの対応付けを記述します。

.. code-block:: python

  def create_modeproc(mode, qe):
      if mode in ["scf", "nscf", "relax", "vc-relax", "bands"]:
          modeproc = QEmode_pw(qe)
      else:
          modeproc = QEmode_generic(qe)
      return modeproc
	

modeごとの変換機能は ``QEmode_base`` の派生クラスとしてまとめられています。
このクラスは
``update_namelist()`` で namelist ブロックの更新と、
``update_cards()`` で cards ブロックのデータ生成を行います。
現在は pw.x の scf、nscf、relax、vc-relax、bands に対応する ``QEmode_pw`` クラスと、変換せずそのまま出力する ``QEmode_generic`` クラスが用意されています。

.. code-block:: python

  class QEmode_base:
      def __init__(self, qe):
      def update_namelist(self, content):
      def update_cards(self, content):


namelist については、空欄の値を結晶構造データ等から生成して代入するほか、格子パラメータなど Structure から決まる値や、他のパラメータとの整合性をとる必要のある値を強制的にセットする場合があります。処理内容はモードごとに個別に対応します。

cards ブロックについては、card の種類ごとに関数を用意し、card名と関数の対応付けを ``card_table`` 変数に列挙します。
基底クラスの ``update_cards()`` では、card名から対応する関数を取得して実行し、card の情報を更新します。もちろん、全く独自に ``update_cards()`` 関数を作成することもできます。

.. code-block:: python

    self.card_table = {
        'CELL_PARAMETERS': generate_cell_parameters,
        'ATOMIC_SPECIES': generate_atomic_species,
        'ATOMIC_POSITIONS': generate_atomic_positions,
        'K_POINTS': generate_k_points,
    }

cardごとの関数は ``src/cif2x/qe/cards.py`` にまとめられており、関数名は ``generate_{card名}`` としています。この関数は card ブロックのパラメータを引数に取り、card名、option、dataフィールドからなる辞書データを返します。


対応する計算モード
----------------------------------------------------------------

``QEmode_pw`` は ``scf``、``nscf``、``relax``、``vc-relax``、``bands`` について
構造カード(``CELL_PARAMETERS``、``ATOMIC_SPECIES``、``ATOMIC_POSITIONS``、
``K_POINTS``)を生成し、``calculation`` をタスクの ``mode`` から設定します。
``&ions``(``relax``/``vc-relax``)・``&cell``(``vc-relax``)ネームリストは
テンプレート/``content`` の記述をそのまま使用します。これらが無いと ``pw.x``
はエラーになるため、必ず記述してください。

``bands`` では ``content`` の K_POINTS オプションに ``crystal_b`` を指定します。
経路は結晶対称性から ``pymatgen`` の高対称 k 経路を用いて生成されます。
他の pw.x タスクと同様に、 ``&control``, ``&system``, ``&electrons`` の各ネームリストは
テンプレートまたは ``content.namelist`` で与え、生成したい構造カードは ``content`` に
列挙する必要があります。例:

.. code-block:: yaml

    structure:
      use_primitive: true
      use_ibrav: false

    tasks:
      - mode: bands
        output_file: bands.in
        template: pw.in_tmpl      # &control, &system, &electrons ネームリストを含む
        content:
          namelist:
            system:
              nbnd: 16        # 伝導帯が必要な場合は明示的に指定
          CELL_PARAMETERS:
          ATOMIC_SPECIES:
          ATOMIC_POSITIONS:
            option: crystal
          K_POINTS:
            option: crystal_b
            line_npoints: 30   # 高対称点間に生成する点数

``line_npoints`` の既定値は 20 です。高対称経路に不連続(分断された区間)が
含まれる場合、各区間終端の整数を 0 にすることで ``bands.x`` が区間をまたいで
補間しないようにします。独自経路を使う場合は ``path`` にラベル列のリストを
与えます(ラベルは自動生成された k 経路に含まれる必要があります)。

Quantum ESPRESSO は既定で占有(価電子)バンドのみを計算するため、伝導帯を
バンド構造に含めるには ``content.namelist.system`` で ``nbnd`` を明示的に指定してください。

高対称点の座標は**標準プリミティブセル**で定義されるため、``bands`` では構造が
そのセルである必要があります。``use_primitive: true``(かつ ``use_ibrav: false``)
を指定してください。慣用セルやスーパーセルなど非標準の入力の場合、``cif2x`` は
警告(「band path: ... the path may be incorrect」)を出力し、サンプリングされる
経路は出力される ``CELL_PARAMETERS`` と一致しません。

pw.x に ``calculation='dos'`` はありません。状態密度計算は
``scf`` -> ``nscf``(密なメッシュ)-> ``dos``(``dos.x`` 入力)の連鎖です。
``dos`` タスクは ``QEmode_generic`` でそのまま出力されるため、``&DOS``
ネームリストはテンプレート/``content`` に記述します:

.. code-block:: yaml

    tasks:
      - mode: scf
        output_file: scf.in
        template: pw.in_tmpl      # &control, &system, &electrons ネームリストを含む
        content:
          CELL_PARAMETERS:
          ATOMIC_SPECIES:
          ATOMIC_POSITIONS:
            option: crystal
          K_POINTS:
            option: automatic
            grid: [8, 8, 8]
      - mode: nscf
        output_file: nscf.in
        template: pw.in_tmpl
        content:
          namelist:
            system:
              occupations: tetrahedra
          CELL_PARAMETERS:
          ATOMIC_SPECIES:
          ATOMIC_POSITIONS:
            option: crystal
          K_POINTS:
            option: automatic
            grid: [16, 16, 16]
      - mode: dos
        template: dos.in_tmpl     # &DOS ネームリストを含む
        output_file: dos.in

生成した入力は ``scf`` → ``nscf`` → ``dos.x`` の順に実行し、各ステップが前の
ステップの保存データを読めるよう、3つすべてで同じ ``prefix`` と ``outdir`` を
指定してください。


トラブルシューティング
----------------------------------------------------------------

cif2x 実行時に起こりやすいエラーと対処法を以下にまとめます。

cif2x は入力ファイルを生成する前に、ターゲット名と ``input.yaml`` を検証します。ターゲット名や入力が不正な場合は、エラーメッセージを1行出力して終了コード ``1`` で終了します (Python のトレースバックは表示されません)。

- **unsupported target '...'**

  ``-t``/``--target`` の値が ``quantum_espresso`` (``qe``, ``espresso``)、 ``vasp``、 ``openmx``、 ``akaikkr``、 ``respack`` のいずれでもない場合に出力されます。メッセージに有効な選択肢が列挙されます。

- **task N: 'mode' is required for target '...'**

  Quantum ESPRESSO または RESPACK 向けの実行で、 ``tasks`` の各要素に ``mode`` が指定されていない場合に出力されます。各 task に ``mode`` (``scf``, ``nscf`` など) を指定してください。

- **task N: 'output_file' is required for target '...'**

  出力ファイル名が必要な task に ``output_file`` が指定されていない場合に出力されます。各 task に出力ファイル名を ``output_file`` で指定してください (``quantum_espresso``、 ``openmx``、 ``akaikkr``、 ``respack`` で必須です)。

- **task N: unknown key '...' for target '...'**

  task に未知のキー (例: ``output_fil`` のようなタイプミス) が含まれる場合に出力されます。メッセージにそのターゲットで許可されるキーが列挙されます。なお、自由記述ブロック (``content``、 ``optional``、 ``structure``) のキーは検証対象外です。

- **pp_file / cutoff_file / pseudo_dir not specified または not found の警告**

  ``optional`` セクションでこれらのパラメータが指定されていない、または指定したパスが存在しない場合に警告が出力されます (処理は継続します)。擬ポテンシャルのインデックスファイル (``pp_file``)、カットオフのインデックスファイル (``cutoff_file``)、擬ポテンシャル格納ディレクトリ (``pseudo_dir``) のパスが正しいか確認してください。 ``pp_file`` は ``ATOMIC_SPECIES`` を生成する場合や ``nbnd`` を空欄にした場合に必須で、未指定のままだと **ATOMIC_SPECIES: optional.pp_file ... is required** または **nbnd is left blank but ...** のエラーで停止します。

- **cutoff information (ecutwfc) not found for element '...'**

  ``ecutwfc`` を空欄 (または未記載) にしたが、その元素について解決できなかった場合に出力されます。 ``cutoff_file`` にその擬ポテンシャルファイルの行がなく、 ``pseudo_dir`` 内の ``.UPF`` ファイルにもカットオフ情報がない状態です。 ``cutoff_file`` に行を追加する、 ``.UPF`` ファイルを配置する、または ``content.namelist.system.ecutwfc`` を明示的に指定してください。 ``ecutrho`` が解決できない場合はエラーにならず、キーを省いて ``pw.x`` のデフォルト (``4 * ecutwfc``) に任せます。

- **cannot determine the valence charge of '...'**

  結晶構造データが磁気モーメントを持ち ``starting_magnetization(i)`` を生成する際に、擬ポテンシャルの価電子数が分からない場合に出力されます。 ``pp_file`` に ``zval`` 列を追加する、 ``pseudo_dir`` に ``.UPF`` ファイルを配置する、または ``content.namelist.system`` で ``starting_magnetization(i)`` を明示的に指定してください。

- **POTCAR: ...** (VASP)

  POTCAR を組み立てられなかった場合に出力されます。 ``pp_file`` に元素の行がない、POTCAR ファイルが読めない、 ``pp_file`` も ``PMG_VASP_PSP_DIR`` も設定されていない、のいずれかです。 ``optional`` セクションを修正するか (ファイルフォーマットの章を参照)、意図的に POTCAR を省く場合は ``content.potcar.skip: true`` を指定してください。

- **the structure has partially occupied (disordered) sites**

  Quantum ESPRESSO、VASP、OpenMX の入力生成には規則構造が必要です。合金は規則配置のスーパーセルモデルにするか、コヒーレントポテンシャル近似で不規則サイトを扱える ``akaikkr`` ターゲットを使用してください。
