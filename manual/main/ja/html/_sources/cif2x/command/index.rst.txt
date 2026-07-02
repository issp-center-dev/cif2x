コマンドリファレンス
================================================================

cif2x
----------------------------------------------------------------

  第一原理計算のための入力ファイルを生成する

書式:

  .. code-block:: bash

    cif2x [-v][-q][--dry-run] -t target input_yaml material.cif
    cif2x [-v][-q][--dry-run] -t target --mp-id ID [--symprec PREC][--api-key-file FILE] input_yaml
    cif2x -h
    cif2x --version

説明:

  input_yaml に指定した入力パラメータファイルと material.cif に指定した結晶構造データを読み込み、第一原理計算プログラム用の入力ファイルを生成します。現在は Quantum ESPRESSO, VASP, OpenMX, AkaiKKR に対応しています。
  以下のオプションを受け付けます。

  - ``-v``

    実行時に表示されるメッセージを冗長にします。複数回指定すると冗長度が上がります。
    
  - ``-q``

    実行時に表示されるメッセージの冗長度を下げます。 ``-v`` の効果を打ち消します。複数回の指定が可能です。

  - ``-t`` *target*

    対象となる第一原理計算プログラムを指定します。 *target* として指定可能な値は以下のとおりです。大文字小文字は区別しません。

    - ``QE``, ``espresso``, ``quantum_espresso``: Quantum ESPRESSO向け入力ファイルを生成します。

    - ``VASP``: VASP向け入力ファイルを生成します。

    - ``OpenMX``: OpenMX向け入力ファイルを生成します。

    - ``AkaiKKR``: AkaiKKR向け入力ファイルを生成します。

    - ``respack``: RESPACK ワークフロー一式を生成します。Quantum ESPRESSO 入力(``scf``/``nscf`` タスク、および任意で ``bands`` タスク)と RESPACK 制御ファイル ``input.in``(``mode: respack`` タスク)を出力します。構造は pymatgen の高対称 k 経路が前提とする**標準プリミティブセル**である必要があります(``structure.use_primitive: true``、``use_ibrav: false`` を指定し、標準プリミティブ構造を与えてください)。非立方系で不一致のセルを与えると、書き出される k 経路が QE のセルと不整合になります(cif2x は警告を出力します)。``N_wannier`` とエネルギー窓は ``content``/テンプレートで与える物理量で、初期ゲスは SCDM(``N_initial_guess = 0``、``respack-wannier-py`` を対象)です。``nscf`` タスクには ``qe2respack`` 用に ``nosym = .true.``/``noinv = .true.`` を設定してください。実行順序: ``scf`` → ``nscf`` → ``qe2respack`` → ``respack``。

  - ``--dry-run``

    生成される入力ファイルをディスクに書き込まず、標準出力に表示します。ファイルやディレクトリを作成せずに生成結果を確認したい場合に便利です。

  - ``input_yaml``

    入力パラメータファイルを指定します。形式は YAML format です。

  - ``material.cif``

    結晶構造データファイルを指定します。形式は CIF の他、pymatgen で扱える形式のファイルを指定できます。``--mp-id`` を指定する場合は省略可能です。

  - ``--mp-id`` *ID*

    ``material.cif`` を読む代わりに、Materials Project の物質 ID *ID* (例: ``mp-149``)の結晶構造を直接取得します。``material.cif`` と ``--mp-id`` はいずれか一方のみを指定してください。API キーは ``--api-key-file``(既定 ``materials_project.key``)、続いて環境変数 ``MP_API_KEY``、pymatgen 設定ファイル(``~/.config/.pmgrc.yaml``)の ``PMG_MAPI_KEY`` の順に、``getcif`` と同じ方法で解決されます。

  - ``--symprec`` *PREC*

    取得した構造を書き出す際の対称性の許容誤差(symprec、既定 ``0.1``、Materials Project に準拠)。``--symprec 0`` で対称性の精密化を無効化します。``--mp-id`` 指定時のみ有効です。

  - ``--api-key-file`` *FILE*

    Materials Project の API キーを格納したファイル(``#`` 以外の各行に1つ、既定 ``materials_project.key``)。``--mp-id`` 指定時のみ有効です。

  .. note::

     ``--mp-id`` は、取得した構造を一時 CIF に書き出して読み直すことで「``getcif`` のあと ``cif2x``」の流れを再現します。Materials Project の最終構造をそのまま使用します(標準セル〔慣用セル〕への変換はしません)。CIF を経由するため、磁気モーメントなどのサイト特性は引き継がれず、無秩序(部分占有)構造は Quantum ESPRESSO 生成で拒否されます。これらは2段階の流れと同じ制限です。

  - ``-h``

    ヘルプを表示します。

  - ``--version``

    バージョン情報を表示します。

