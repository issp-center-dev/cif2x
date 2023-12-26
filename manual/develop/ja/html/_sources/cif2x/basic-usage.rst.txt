インストールと基本的な使い方
================================================================

**必要なライブラリ・環境**

  HTP-tools に含まれる第一原理計算入力ファイル生成ツール cif2x を利用するには、以下のプログラムとライブラリが必要です。

  - python 3.x
  - pymatgen モジュール
  - ruamel.yaml モジュール
  - f90nml モジュール
  - qe-tools モジュール
  - numpy モジュール
  - pandas モジュール
  - monty モジュール

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

  実行プログラム ``cif2x`` がインストールされます。

**ディレクトリ構成**

  ::

     .
     |-- LICENSE
     |-- README.md
     |-- pyproject.toml
     |-- docs/
     |   |-- ja/
     |   |-- tutorial/
     |-- src/
     |   |-- cif2x/
     |       |-- __init__.py
     |       |-- main.py
     |       |-- cif2struct.py
     |       |-- struct2qe.py
     |       |-- qe/
     |       |   |-- __init__.py
     |	     |   |-- calc_mode.py
     |	     |   |-- cards.py
     |	     |   |-- content.py
     |	     |   |-- qeutils.py
     |	     |   |-- tools.py
     |       |-- struct2vasp.py
     |       |-- struct2openmx.py
     |       |-- openmx/
     |       |   |-- __init__.py
     |       |   |-- vps_table.py
     |       |-- utils.py
     |-- sample/


**基本的な使用方法**

  cif2xは第一原理計算プログラムのための入力ファイルを生成するツールです。入力パラメータを雛形として、物質の種類や計算条件によって変わる箇所を結晶構造データなどから構成します。現在は Quantum ESPRESSO, VASP, および OpenMX の入力ファイル形式に対応しています。

  #. 入力パラメータファイルの作成

      cif2xを使用するには、まず、生成する入力ファイルの内容を記述したパラメータファイルをYAML形式で作成します。詳細についてはファイルフォーマットの章を参照してください。

  #. 結晶構造ファイルと擬ポテンシャルファイルの配置

      対象となる物質の結晶構造を記述したファイルを用意します。ファイル形式は CIF または pymatgen が扱える POSCAR や xfs 形式に対応しています。

      Quantum ESPRESSO の場合、利用する擬ポテンシャルファイルと、CSV形式のインデックスファイルを配置します。擬ポテンシャルファイルの配置先などは入力パラメータファイル内に指定します。

      VASP の場合、擬ポテンシャルファイルの格納場所を ``~/.config/.pmgrc.yaml`` ファイルに記述するか環境変数にセットします。入力パラメータファイル内で指定することもできます。

  #. コマンドの実行

      作成した入力パラメータファイルおよび結晶構造データファイルを入力としてcif2xプログラムを実行します。Quantum ESPRESSO用の入力ファイルを生成する場合はターゲットオプションに ``-t QE`` を指定します。VASPの場合は ``-t VASP``, OpenMX の場合は ``-t OpenMX`` を指定します。

      .. code-block:: bash

          $ cif2x -t QE input.yaml material.cif

