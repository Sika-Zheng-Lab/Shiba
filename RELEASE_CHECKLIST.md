# v1.0.0 リリース前チェック

確認日: 2026-11-13 (America/Los_Angeles)

**判定: 本番公開前に下記の未完了項目への対応が必要。** developへの統合と検証を実施した。mainへのマージ、v1.0.0タグ作成、GitHub Release公開、Docker Hubへの配布は行っていない。

## 統合と今回の修正

- [x] `feature/event-centric-betabinomial` の6コミットをdevelopへ競合なくマージし、origin/developへpush。マージコミット: `81ed5e05f57da4892fd82b9de5a51c13127df94d`。
- [x] CIで判明したPython標準の`test`パッケージとの名前衝突を`test/__init__.py`で解消。
- [x] Arialがない環境ではDejaVu Sansを許容する仕様に合わせてSVG出力テストを修正し、フォールバックのSVGも検証。
- [x] CIに不足していた`pysam` / `pyyaml`を追加し、Rバックエンドを使うジョブにYAMLランチャーの統合テストを追加。
- [x] 小規模GTFで一部のイベント種類が0件になるとbulk Fisher解析が`KeyError: event_id`で停止する不具合を修正。イベント選択マスクを明示的なbool型にし、全8種類の空テーブルで列が維持されることを回帰テストで確認。

コード修正のコミットは`0b1ca60`（CI）と`ccb2f03`（空イベント）。この記録の検証対象コードは`ccb2f03`。

## 実行結果

| 確認 | 結果・範囲 |
| --- | --- |
| Pythonテスト一式 | 243件、失敗0、Rバックエンド未導入によるskip 13件。残るR検証は下段に別記 |
| R/glmmTMB実バックエンド | 最終コード`ccb2f03`のGitHub Actionsで35件すべて成功、skipなし。ローカルの13件のskipを合格扱いにせず、実バックエンドで別途検証した |
| GitHub Actions | 最終コード`ccb2f03`の[run](https://github.com/Sika-Zheng-Lab/Shiba/actions/runs/43491174320)で通常テスト・R統合とも成功。Codecovへのアップロードも成功 |
| `mkdocs build --strict` | 成功。ローカルリンク・アンカーを含むMkDocs設定上の検証 |
| Python構文 / CLI | 追跡対象Pythonファイルのコンパイル、`shiba.py` / `scshiba.py` / `src/psi.py` / `src/scpsi.py`の`--help`成功 |
| Python依存関係 | `pip check`成功 |
| SnakeShiba dry-run | Fisher 16ジョブ、factorial 57ジョブのDAG構築成功。外部コマンドの実行成功を意味しない |
| scShiba + Tosa 1.0.0 | 合成BAMの統合テスト成功。UMI重複除去、ライブラリ横断のグループ集計、RI PSI、全8種類の出力を確認 |
| SnakeScShiba + Tosa 1.0.0 | 合成BAMで全8ジョブ成功。期待junction counts、HTML、report.jsonを確認 |
| bulk Fisher + Tosa 1.0.0 | 合成BAM、既知GTF、`--start-step 2 --mame -p 2`で成功。空イベントの修正後、全8種類のPSI出力まで完走 |
| Fisher Excel出力 | **失敗**。空のSE結果で`styleframe==4.2`が`IndexError: row: 2 is out of rows range`を送出。TSV生成後にExcelで停止する |
| Docker / Apptainer | この環境にランタイムがなく、イメージのbuild/runは未検証 |
| bulkの全工程 | StringTieアセンブリ、featureCounts、DESeq2を含む全工程・実データ検証は今回未実施 |

ローカル環境はPython 3.12.15、pandas 2.3.3、NumPy 2.5.4、Snakemake 9.27.1、Tosa 1.0.0。生ログ、fixture、依存パッケージ一覧は`/tmp/shiba-v1-release-check/`に保存した。一時ディレクトリなので、長期保存が必要な場合は移すこと。

R統合CIはUbuntuの`glmmTMB 1.1.14+dfsg-2`を使用。テストは成功したが、glmmTMBのビルド時TMB 1.9.17と実行時1.9.19の不一致警告が出た。配布用環境では整合した組合せを固定して検証する。ローカルRの追加ビルドは、実バックエンドCIの成功確認後に中止した。

Excelの失敗ログは`bulk-excel.log`。今回生成したfixtureを使う再現コマンド:

```bash
.venv/bin/python src/psi.py \
  /tmp/shiba-v1-release-check/bulk-smoke-output/junctions/junctions.bed \
  /tmp/shiba-v1-release-check/bulk-smoke-output/events \
  /tmp/shiba-v1-release-check/bulk-excel \
  -g /tmp/shiba-v1-release-check/bulk-experiment.tsv \
  -r Ref -a Alt -m 1 -t False --excel True -p 2
```

異常終了後に`results.xlsx`が存在しても成功とは判定できない。今回のファイルにはSEシートしかなく、正常な全8シート出力にはなっていない。

## 公開前に必須の作業

| 優先度 | 作業 | 根拠・完了条件 |
| --- | --- | --- |
| P0 | VERSIONとリリース日を確定する | `VERSION`は`v0.8.2`、CHANGELOGの日付は`2026-XX-XX`。VERSIONを`v1.0.0`に更新し、CLI・README・配布物・設定例を一致させる。現状でmainへマージすると既存v0.8.2タグの作成を試みる |
| P0 | 空イベントのExcel出力を修正する | `src/lib/shibalib.py`の`save_excel` / `save_excel_sc`は空シートでも`B2`をfreezeする。bulkの`--excel True`でSEが0件だと実際に停止した。ヘッダーのみの全8シートを保存できることを検証し、他シートの例外を握りつぶす処理も見直す |
| P0 | バージョン更新の一括置換を修正する | `.github/workflows/bump-version.yaml`と`update-version.yaml`は文中の全バージョンを置換するため、「v0.8.2には新バックエンドがない」が「v1.0.0にはない」に変わる。旧版に言及する移行説明を保ち、更新対象を明示する。factorial用設定例も更新対象に含める |
| P0 | 公開処理の順序と再実行を整える | `release-docker-build-push.yaml`はタグ・Releaseを公開した後でDockerをbuild/pushする。イメージとテストの成功を公開条件にし、既存タグ/Releaseがある場合の再実行を扱う。対象SHA、VERSION形式、タグ未使用も事前検証する |
| P1 | ドキュメント公開の条件を修正する | `publish-mkdocs.yaml`はmain向けPRが閉じると未マージでも動く。`merged == true`の条件を付け、厳密ビルドも公開前のCIで確認する |
| P1 | CHANGELOGと導入案内の矛盾を解消する | CHANGELOGには削除済み`beta_binomial: True`の追加説明と、TosaのBioconda配布前提と食い違うCargo案内が残る。READMEの`docker pull shiba:v0.8.2`は配布先`naotokubota/shiba`とも異なる。development限定の説明、source checkoutの選択手順もv1.0.0向けに見直す |
| P1 | 配布用Dockerイメージを実際にbuild/runする | 両DockerfileでPython/R/Tosa/Excel、シンボリックリンク経由のVERSION解決、bulk/scのsmoke testを確認。`.dockerignore`がなく`COPY . /opt/Shiba`するため、`.git`、`.venv`、検証出力等を除外する。配布対象CPUアーキテクチャと検証した依存バージョンを記録する |
| P1 | Conda配布を準備する | ShibaとMameShibaのBiocondaレシピについて、v1.0.0ソース、checksum、Tosa、R/glmmTMB/jsonlite、新規Rファイル・テンプレート・VERSIONの同梱を確認。レシピはこのリポジトリ内になく、今回は変更・検証していない |
| P1 | 配布候補の統合検証を完了する | 下の検証チェックリストを満たし、最終SHAのCIをgreenにする。Docker Hub資格情報と公開権限の有効性も実際の配布準備で確認する |

P0は現状のまま公開操作を進められない項目、P1は正式版の配布前に完了させる項目。

## 残る検証と移行案内

- [ ] 配布予定の依存環境で、bulk Fisher / factorial、scShiba、SnakeShiba、SnakeScShibaをそれぞれ完走する。今回の小規模smoke testとbulkのdry-runを全工程検証の代わりにしない。
- [ ] RバックエンドとTosaの必須テストがCIでskipにならないようにする。現行CIにはTosa導入ジョブがなく、single-cell BAM統合テストはskipになる。
- [ ] 代表的な実データで、全8イベント、minus-strand MXE、RI、短鎖/長鎖RNA-seq、BAM/CRAM、single-cellのCB/UBタグとbarcode groupを確認する。
- [ ] beta-binomialについて、係数/contrast、交互作用、予測PSI、フィルタとFDRの対象、収束失敗時の扱いを確認する。大規模データの時間・ピークメモリと既知の制約を記録する。
- [ ] HTMLをブラウザで開き、オフライン表示・Plotly操作・PNG/PDF/SVG・Excelの内容を確認する。自動テストは生成ファイルと埋め込み内容を検証するが、ブラウザの操作確認は今回未実施。
- [ ] リリースノートで、Tosaへの移行に伴う再集計、minus-strand MXEの定義変更、single-cellの新experiment table、旧`beta_binomial`キー/CLIの廃止、R依存、factorialはbulk向けであることを明示する。
- [ ] VERSION更新ワークフローを整理する。`update-version.yaml`には旧`::set-output`、既に動的読込へ移行したVERSION文字列の置換、ステップだけを終了して後続ジョブ処理を止めないmerge判定が残る。

## リリース実行の順序

1. 上記P0/P1対応をdevelopへ反映し、実際のリリース日とv1.0.0の表記を確定する。
2. 同じ候補SHAでCI、厳密ドキュメントビルド、配布イメージの統合テストを成功させる。
3. develop → mainのPRに確定した変更点・破壊的変更・検証結果を記載する。現在の自動処理では、このPRのマージが本番公開のトリガーになる。
4. 公開後にv1.0.0タグ/Release、Dockerの`v1.0.0`/`latest`、公開ドキュメントの整合を確認する。
5. BiocondaのShiba/MameShiba配布を確認し、新規環境でインストールと最小実行を再確認する。
