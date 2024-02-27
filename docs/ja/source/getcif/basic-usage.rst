インストールと基本的な使い方
================================================================

**必要なライブラリ・環境**

  HTP-tools に含まれるCIFデータ取得ツール getcif を利用するには、以下のプログラムとライブラリが必要です。

  - python 3.x
  - pymatgen モジュール
  - ruamel.yaml モジュール
  - mp-api モジュール

**ソースコード配布サイト**

  - GitHubリポジトリ `https://github.com/issp-center-dev/getcif <https://github.com/issp-center-dev/getcif>`_

**ダウンロード方法**

  gitを利用できる場合は、以下のコマンドでgetcifをダウンロードできます。

  .. code-block:: bash

    $ git clone https://github.com/issp-center-dev/getcif.git

**インストール方法**

  getcifをダウンロード後、以下のコマンドを実行してインストールします。getcifが利用するライブラリも必要に応じてインストールされます。

  .. code-block:: bash

     $ cd ./getcif
     $ python3 -m pip install .

  実行プログラム ``getcif`` がインストールされます。

**ディレクトリ構成**

  ::

     .
     |-- LICENSE
     |-- README.md
     |-- pyproject.toml
     |-- docs/
     |   |-- ja/
     |   |-- en/
     |   |-- tutorial/
     |-- src/
     |   |-- getcif/
     |       |-- __init__.py
     |       |-- main.py
     |-- sample/


**基本的な使用方法**

  getcifは物質材料データベースから結晶構造データ等を取得するツールです。現在は Materials Project からのデータ取得に対応しています。物質の組成や対称性、バンドギャップなどの物性値をもとにデータベースを検索し、データを取得することができます。

  #. APIキーの取得

      Materials Project のデータを利用するには `Materials Project のウェブサイト <https://next-gen.materialsproject.org/>`_ でユーザ登録が必要です。登録すると、プログラム等からデータを検索するための APIキーを取得できます。APIキー(文字列)は、pymatgen のツールを用いて設定ファイルにセットするか、環境変数に設定します。あるいは、入力パラメータファイルで指定したファイルに APIキーを書き込みます。

  #. 入力パラメータファイルの作成

      getcifを使用するには、検索条件と取得するデータの内容を記述した YAML形式のパラメータファイルを作成します。詳細についてはファイルフォーマットの章を参照してください。

  #. コマンドの実行

      作成した入力パラメータファイルを入力としてgetcifプログラムを実行します。

      .. code-block:: bash

          $ getcif input.yaml

      データベースに接続してデータを取得します。取得したデータは、material ID をディレクトリ名とするディレクトリ内に出力されます。結晶構造データは CIF 形式で保存されます。
