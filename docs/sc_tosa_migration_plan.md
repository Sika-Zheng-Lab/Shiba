# scShiba / SnakeScShiba の Tosa 移行計画

## 現状と目標

- 現行の実験表は `barcode`（`barcode`, `group` 列を持つ TSV）と `SJ`（STARsolo `Solo.out/SJ/raw`）を参照する。`src/sc2junc.py` は STARsolo の疎行列を読み、バーコードごとの junction 数を group に合計し、`chr`, `start`, `end`, `ID`, group 列から成る `junctions.bed` を作る。
- `src/scpsi.py` と `snakescshiba.smk` は SE、FIVE、THREE、MXE、MSE、AFE、ALE の 7 種類だけを処理する。`gtf2event.py` は `EVENT_RI.txt` も作るが、単一細胞ワークフローは現在これを使わない。
- Tosa single-cell モードで各ライブラリの BAM/CRAM から junction と exon–intron 境界を数え、同じ group 列の `junctions.bed` に統合する。RI を 8 番目の解析対象に加える。既存 7 種類の出力形式と group の意味は維持する。

## Tosa 1.0.0 の入出力契約

```bash
tosa single -c selected_barcodes.tsv -g annotation.gtf \
  -a 8 -b 1 -m 20 -M 500000 -p 8 input.bam output_prefix
```

`input.bam` は CRAM も指定できる。`-c` はヘッダーなし 1 列のバーコード一覧で、Shiba の `barcode,group` 表をそのまま渡せない。Tosa は BAM/CRAM の `CB` タグで細胞を識別し、`UB` タグがあれば同一 junction/境界上の同一 `(CB, UB)` を重複排除する。`UB` がない場合は read 単位のカウントとなるため、比較前にタグの存在率を調べる。`-g` はイベント生成に使う GTF を渡す。現行 scShiba は STARsolo junction の strand を統合しているため、最初の比較では Tosa の `-s` を省略し、設定で `XS` / `RF` / `FR` を選べるようにする。

Tosa single-cell の junction 出力は `{prefix}_matrix.mtx.gz`（行: feature、列: barcode）、`_barcodes.tsv.gz`、`_features.tsv.gz`（座標と strand）、`_junction_barcodes.tsv.gz`（`Feature`, `Strand`, `Barcode`, `Count`）。GTF 指定時は `{prefix}_boundary_matrix.mtx.gz`、`_boundary_barcodes.tsv.gz`、`_boundary_features.tsv.gz`（座標、`5p`/`3p`、strand）、`_boundary_barcodes_detail.tsv.gz`（`Boundary`, `Type`, `Strand`, `Barcode`, `Count`）も作る。詳細 TSV を gzip ストリームとして読み、行列全体をメモリに展開せず group 別に合計する。Snakemake では使う出力を明示して欠落を検知する。

Tosa junction の `chr:S-E` は intron 内部の 1-based inclusive 座標であり、現行 scShiba の ID は `chr:(S-1)-(E+1)`。Tosa 境界は 0-based half-open `chr:start-end` で、RI 用 ID は `5p → chr:(start+b)-(start+b+1)`、`3p → chr:(end-b)-(end-b+1)` とする（`b` は `-b` の設定値）。`EVENT_RI.txt` に含まれる境界だけを統合し、検出されない境界もゼロ行として出力する。染色体名はイベント GTF と照合し、既存 scShiba が持つ contig 名を不用意に書き換えない。

## 実装手順

1. **入力仕様と検証**: 単一細胞実験表を `sample`, `alignment`, `barcode` の 3 列に変更する。`sample` はライブラリごとの一意な名前、`alignment` は BAM/CRAM、`barcode` は現行の `barcode,group` TSV。`src/lib/general.py` の sc 検証を更新し、ファイル存在、重複 sample・barcode/group の矛盾、参照 group、Tosa の anchor/intron/strand 設定を確認する。旧 `SJ` のみの表は移行方法を示すエラーにする。
2. **Tosa 呼び出しと集計**: `src/sc2junc.py` を Tosa single-cell の実行器・出力変換器に更新する。各 barcode/group TSV から重複のない 1 列 whitelist を生成し、各 `alignment` に対して `-c`, `-g`, `-a`, `-b`, `-m`, `-M`, `-p` と必要なら `-s` を渡す。詳細 TSV のバーコードをそのライブラリの group 表へ結合し、junction と RI 境界を group ごとに合計する。同名バーコードが別ライブラリにあっても各行の group 表を個別に適用する。`junctions.bed` の既存列とゼロ埋めを維持する。
3. **両ワークフロー**: `scshiba.py` は Step 2 に GTF、RI イベント、Tosa 設定とスレッド数を渡す。`snakescshiba.smk` は `sample` ごとの Tosa rule と統合 rule に分け、BAM/CRAM、barcode TSV、GTF、`EVENT_RI.txt` を依存関係にする。Tosa の junction/境界出力を明示し、両経路が同じ変換関数を使う。
4. **RI 解析**: `src/lib/shibalib.py` の既存 `event_for_analysis_ri`, `ri`, `ri_ind`, `diff_ri` を `src/scpsi.py` に接続する。`read_events_sc`、イベント定義、`PSI_RI.txt`、summary、`onlypsi` の PSI 行列、`save_excel_sc` に RI を加える。`snakescshiba.smk` の `rule all` / `gtf2event` / `scpsi` の入出力も 8 種類に更新する。
5. **利用者向け資料**: 単一細胞の実験表・設定例、BAM/CRAM の `CB` / `UB` 要件、Tosa 依存、RI の追加、既存結果との比較上の注意を README・マニュアル・CHANGELOG に反映する。

## 検証と完了条件

1. Tosa 同梱の synthetic BAM/GTF/barcodes を `single -b 1` で実行し、junction（2 feature × 3 barcode）と境界（4 feature × 3 barcode）の行列・詳細 TSV の対応、座標変換、group 集計、RI のゼロ行を確認する。実測では `AAAA-1` / `BBBB-1` / `CCCC-1` の 3 barcode、junction 2 feature、境界 4 feature が出力された。
2. 複数ライブラリ、同名バーコード、同じ座標の ± strand、未検出 RI 境界、バーコード表にない CB、空の結果、CRAM、タグ欠損を単体・小規模統合テストに含める。全 8 種類の PSI と `onlypsi`・差次的解析・Excel 出力を検証する。
3. 旧 `test_scshiba.bash` / `test_SnakeScShiba.bash`、`test_config_scShiba.yaml`、`test_experiment_table_scShiba.tsv` と保存済み結果を比較基準にする。対応する元の BAM/CRAM が必要。現在の `/rhome/naotok/bigdata/Shiba/test_scshiba/` には STARsolo の SJ 行列と barcode 表だけが見つかり、BAM/CRAM は確認できていない。元データが得られたら新形式の実験表を作り、共通 junction ID の group カウントと 7 種類の既存 PSI・差次的判定の変化を集計する。Tosa と STARsolo のフィルタ・UMI 集計差は記録し、完全一致を前提にしない。
4. RI は旧 scShiba に比較値がないため、イベントごとの `intron_a` と 2 境界のカウント、RI PSI の式、最小 read 閾値、group 別差次的判定を手計算可能な例と bulk RI 関数の結果で検証する。scShiba と SnakeScShiba で同一入力から同じ `junctions.bed` と 8 種類の PSI を生成し、SnakeScShiba の dry-run と実行を行う。

完了条件は、STARsolo `SJ` がなくても BAM/CRAM と barcode/group 表から両ワークフローが終了し、8 種類のイベントに結果を出すこと、既存 7 種類の結果との差分と Tosa 由来の理由を記録すること、必要な Tosa 版と設定・入力を再現できることである。

## 合成データ検証（2026-09-28）

実データの BAM/CRAM がないため、`test/make_sc_tosa_fixture.py` で 2 ライブラリの座標ソート済み BAM と index、GTF、barcode/group 表、実験表、設定ファイルを生成した。同じ `AAAA-1` バーコードを 2 ライブラリで異なる group に割り当て、`UB` 重複、対象外の `CB`、spliced read と RI 両境界を含めた。生成物は `/rhome/naotok/bigdata/tmp/Shiba_sc_tosa_fixture/` に保存した。

scShiba の全 3 ステップ、SnakeScShiba の全 7 ジョブ、`scpsi.py --onlypsi True` が完了した。両ワークフローの `junctions.bed` と 8 種類の `PSI_*.txt` はバイト単位で一致した。`chr1:1200-1500` の junction 数は Ref=2、Alt=3、両 RI 境界はそれぞれ Ref=1、Alt=2。RI PSI は Ref=1/3、Alt=2/5 で、手計算と一致した。Tosa 同梱の synthetic CRAM も `sc2junc.py run` で読み、2 barcode の junction count（6、2）を確認した。回帰テストは `test/test_sc_tosa_integration.py` に追加した。

本データの BAM/CRAM が提供されたら、旧 STARsolo SJ と 7 種類の PSI の比較を再実行する。今回の合成テストでは Excel 出力と本番サイズでのメモリ使用量は未測定である。
