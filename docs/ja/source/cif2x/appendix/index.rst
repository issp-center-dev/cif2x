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
テンプレート/``content`` の記述をそのまま使用します。

``bands`` では ``content`` の K_POINTS オプションに ``crystal_b`` を指定します。
経路は結晶対称性から ``pymatgen`` の高対称 k 経路を用いて生成されます。例:

.. code-block:: yaml

    tasks:
      - mode: bands
        output_file: bands.in
        content:
          system:
            nbnd: 16        # 伝導帯が必要な場合は明示的に指定
          K_POINTS:
            option: crystal_b
            line_npoints: 30   # 高対称点間に生成する点数

``line_npoints`` の既定値は 20 です。高対称経路に不連続(分断された区間)が
含まれる場合、各区間終端の整数を 0 にすることで ``bands.x`` が区間をまたいで
補間しないようにします。独自経路を使う場合は ``path`` にラベル列のリストを
与えます(ラベルは自動生成された k 経路に含まれる必要があります)。なお
``bands`` は当面 ``use_ibrav: false`` を対象とします。

pw.x に ``calculation='dos'`` はありません。状態密度計算は
``scf`` -> ``nscf``(密なメッシュ)-> ``dos``(``dos.x`` 入力)の連鎖です。
``dos`` タスクは ``QEmode_generic`` でそのまま出力されるため、``&DOS``
ネームリストはテンプレート/``content`` に記述します:

.. code-block:: yaml

    tasks:
      - mode: scf
        output_file: scf.in
        content: { K_POINTS: { option: automatic, grid: [8, 8, 8] } }
      - mode: nscf
        output_file: nscf.in
        content: { K_POINTS: { option: automatic, grid: [16, 16, 16] } }
      - mode: dos
        template: dos.in_tmpl     # &DOS ネームリストを含む
        output_file: dos.in


トラブルシューティング
----------------------------------------------------------------

cif2x 実行時に起こりやすいエラーと対処法を以下にまとめます。

cif2x は入力ファイルを生成する前に、ターゲット名と ``input.yaml`` を検証します。ターゲット名や入力が不正な場合は、エラーメッセージを1行出力して終了コード ``1`` で終了します (Python のトレースバックは表示されません)。

- **unsupported target '...'**

  ``-t``/``--target`` の値が ``quantum_espresso`` (``qe``, ``espresso``)、 ``vasp``、 ``openmx``、 ``akaikkr`` のいずれでもない場合に出力されます。メッセージに有効な選択肢が列挙されます。

- **task N: 'mode' is required for target 'quantum_espresso'**

  Quantum ESPRESSO 向けの実行で、 ``tasks`` の各要素に ``mode`` が指定されていない場合に出力されます。各 task に ``mode`` (``scf``, ``nscf`` など) を指定してください。

- **task N: 'output_file' is required for target '...'**

  出力ファイル名が必要な task に ``output_file`` が指定されていない場合に出力されます。各 task に出力ファイル名を ``output_file`` で指定してください (``quantum_espresso``、 ``openmx``、 ``akaikkr`` で必須です)。

- **task N: unknown key '...' for target '...'**

  task に未知のキー (例: ``output_fil`` のようなタイプミス) が含まれる場合に出力されます。メッセージにそのターゲットで許可されるキーが列挙されます。なお、自由記述ブロック (``content``、 ``optional``、 ``structure``) のキーは検証対象外です。

- **pp_file / cutoff_file / pseudo_dir not specified または not found の警告**

  ``optional`` セクションでこれらのパラメータが指定されていない、または指定したパスが存在しない場合に警告が出力されます (処理は継続します)。擬ポテンシャルのインデックスファイル (``pp_file``)、カットオフのインデックスファイル (``cutoff_file``)、擬ポテンシャル格納ディレクトリ (``pseudo_dir``) のパスが正しいか確認してください。なお、 ``cutoff_file`` / ``pseudo_dir`` は警告のみで処理が継続しますが、 ``pp_file`` は ``ATOMIC_SPECIES`` の生成や ``nbnd`` の自動設定で実際に参照されるため事実上必須であり、未指定のままだと処理の後段で ``RuntimeError`` となります。

- **カットオフが 0.0 になる**

  ``ecutwfc`` / ``ecutrho`` を空欄にして自動取得させる場合、(その元素について ``pp_file`` の対応付けが存在する前提で) 対応する ``.UPF`` ファイルが ``pseudo_dir`` に見つからず、かつ ``cutoff_file`` にも該当エントリがないと、カットオフの値が ``0.0`` にフォールバックします。擬ポテンシャルファイルやカットオフのインデックスエントリが揃っているか確認してください。(なお、 ``pp_file`` にその元素の行自体が無い場合は ``0.0`` フォールバックではなく ``KeyError`` で停止します。)
