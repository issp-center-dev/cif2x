インストールと基本的な使い方
================================================================

**必要なライブラリ・環境**

  HTP-tools に含まれる第一原理計算入力ファイル生成ツール cif2x および CIFデータ取得ツール getcif を利用するには、以下のプログラムとライブラリが必要です。
  (これらは ``pyproject.toml`` に記載されており、インストール時に自動的に導入されます。)

  - Python 3.11 以上
  - numpy モジュール
  - pandas モジュール
  - pymatgen モジュール
  - f90nml モジュール
  - beautifulsoup4 モジュール
  - mp_api モジュール
  - ruamel.yaml モジュール

  AkaiKKR 向けの入力ファイルを生成する場合は、上記に加えて AkaiKKRPythonUtil モジュールが必要です。これは自動ではインストールされないため、別途インストールします (下記参照)。

  本バージョンは Python 3.11〜3.14、pymatgen 2026.10、f90nml 1.5、mp-api 0.47 で動作を確認しています。AkaiKKRPythonUtil 2023.2.1 は pymatgen と spglib の非推奨インタフェースに依存しており、pymatgen 2026.10 では動作しますが、これらのライブラリの将来の版で動かなくなる可能性があります。

**ソースコード配布サイト**

  - `GitHubリポジトリ <https://github.com/issp-center-dev/cif2x>`_

**ダウンロード方法**

  gitを利用できる場合は、以下のコマンドでcif2xをダウンロードできます。

  .. code-block:: bash

    $ git clone https://github.com/issp-center-dev/cif2x.git

**インストール方法**

  cif2xをダウンロード後、以下のコマンドを実行してインストールします。cif2xが利用するライブラリも必要に応じてインストールされます。

  .. code-block:: bash

     $ cd ./cif2x
     $ python3 -m pip install .

  実行プログラム ``cif2x``, ``getcif``, ``cif2x-pp-cutoff`` がインストールされます。

  なお、AkaiKKRPythonUtil モジュールは別途インストールが必要です。
  以下の手順で `配布サイト <https://github.com/AkaiKKRteam/AkaiKKRPythonUtil>`_ からソースコードを取得しインストールします。また、必要な seaborn モジュールもインストールしておきます。

  .. code-block:: bash

     $ git clone https://github.com/AkaiKKRteam/AkaiKKRPythonUtil.git
     $ cd ./AkaiKKRPythonUtil/library/PyAkaiKKR
     $ python3 -m pip install .
     $ python3 -m pip install seaborn


**ディレクトリ構成**

  ::

     .
     |-- CHANGELOG.md
     |-- LICENSE
     |-- README.md
     |-- pyproject.toml
     |-- docs/
     |   |-- en/
     |   |-- ja/
     |   |-- tutorial/
     |-- sample/
     |   |-- cif2x/
     |   |-- getcif/
     |-- src/
     |   |-- cif2x/
     |   |   |-- __init__.py
     |   |   |-- main.py
     |   |   |-- input_validator.py
     |   |   |-- cif2struct.py
     |   |   |-- mp_source.py
     |   |   |-- utils.py
     |   |   |-- pp_cutoff.py
     |   |   |-- struct2qe.py
     |   |   |-- qe/
     |   |   |   |-- __init__.py
     |   |   |   |-- calc_mode.py
     |   |   |   |-- cards.py
     |   |   |   |-- content.py
     |   |   |   |-- qeutil.py
     |   |   |   |-- tools.py
     |   |   |-- struct2vasp.py
     |   |   |-- struct2openmx.py
     |   |   |-- openmx/
     |   |   |   |-- vps_table.py
     |   |   |-- struct2akaikkr.py
     |   |   |-- akaikkr/
     |   |   |   |-- make_input.py
     |   |   |   |-- read_input.py
     |   |   |   |-- run_cif2kkr.py
     |   |   |-- struct2respack.py
     |   |-- getcif/
     |       |-- __init__.py
     |       |-- main.py
     |       |-- mp.py
     |-- tests/


**基本的な使用方法**

  cif2xは第一原理計算プログラムのための入力ファイルを生成するツールです。入力パラメータを雛形として、物質の種類や計算条件によって変わる箇所を結晶構造データなどから構成します。現在は Quantum ESPRESSO, VASP, OpenMX, および AkaiKKR の入力ファイル形式に対応しているほか、RESPACK ワークフロー (Quantum ESPRESSO と RESPACK の入力) を生成できます。

  #. 入力パラメータファイルの作成

      cif2xを使用するには、まず、生成する入力ファイルの内容を記述したパラメータファイルをYAML形式で作成します。詳細についてはファイルフォーマットの章を参照してください。

  #. 結晶構造ファイルと擬ポテンシャルファイルの配置

      対象となる物質の結晶構造を記述したファイルを用意します。ファイル形式は CIF または pymatgen が扱える POSCAR や xsf 形式に対応しています。

      Quantum ESPRESSO の場合、利用する擬ポテンシャルファイルと、CSV形式のインデックスファイルを配置します。擬ポテンシャルファイルの配置先などは入力パラメータファイル内に指定します。

      VASP の場合、擬ポテンシャルファイルの格納場所を ``~/.config/.pmgrc.yaml`` ファイルに記述するか環境変数にセットします。入力パラメータファイル内で指定することもできます。

  #. コマンドの実行

      作成した入力パラメータファイルおよび結晶構造データファイルを入力としてcif2xプログラムを実行します。Quantum ESPRESSO用の入力ファイルを生成する場合はターゲットオプションに ``-t QE`` を指定します。VASPの場合は ``-t VASP``, OpenMX の場合は ``-t OpenMX``, AkaiKKR の場合は ``-t akaikkr``, RESPACK ワークフローの場合は ``-t respack`` を指定します。

      .. code-block:: bash

          $ cif2x -t QE input.yaml material.cif

