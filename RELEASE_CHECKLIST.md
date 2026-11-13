# v1.0.0 リリース前チェック

確認日: 2026-11-13 (America/Los_Angeles)

`VERSION`とリリース日の確定は依頼に従い保留する。`VERSION`は`v0.8.2`、CHANGELOGの新バージョンの日付は`2026-XX-XX`のまま。本番mainへのマージ、タグ作成、GitHub Release公開、Docker Hubへのpush、Biocondaへの提出は行っていない。

## 修正状況

| 項目 | 対応 |
| --- | --- |
| VERSIONとリリース日 | **保留（今回の対象外）**。不一致・未確定の日付ではrelease preflightが失敗し、公開処理へ進まない |
| developへの統合 | feature/event-centric-betabinomialを`81ed5e0`で統合。修正はdevelopへ反映 |
| 空イベントのExcel | 全8シートが空でも保存可能。openpyxlに統一し、StyleFrameの追加インストールを不要にした。途中で失敗したファイルは公開せず、既存ファイルを保持する |
| バージョン更新 | `scripts/update_version.py`で更新対象を明示。旧版に関する説明・CHANGELOG・リリース日は置換しない。factorial設定例も対象。旧set-outputとpush時の自動一括置換を廃止 |
| 公開処理 | 同じSHAのテスト・厳密ドキュメントビルド・両Dockerイメージの検証を公開前に必須化。実際に検証した本番イメージをartifact経由で公開する |
| 再実行 | 既存タグ/イメージのSHAが異なれば停止。同じSHAなら再利用。Releaseはdraftへ成果物を登録した後に公開。mainの先端でなくなった古い候補ではlatestを更新しない |
| ドキュメント公開 | mainへのマージ済みPRかつRelease公開成功後のみ実行。未マージで閉じたPRでは公開しない |
| CHANGELOG・導入案内 | 削除済みbeta_binomialキーの説明、Cargo必須の案内、誤ったDocker namespaceを修正。bulk回帰とsingle-cell Fisherの対応範囲を明記 |
| Docker・依存環境 | 共通`environment.yml`、linux/amd64、`.dockerignore`を導入。glmmTMB 1.1.13 / TMB 1.9.19を固定し、起動時の不一致警告を解消。ユーザーsite-packagesの混入を防止 |
| クリーン環境の不具合 | `sample_...`検体名をヘッダーと誤認するbulk処理を修正。PYTHONPATH未設定で止まる両Snakemakeレポート生成を修正 |
| Conda配布準備 | `scripts/prepare_conda_recipes.py`でShiba/MameShibaを生成。検証済みソースarchiveのSHA256を読み込み、Releaseへレシピも添付する。最終版のversion/checksumはVERSION確定後に生成する |
| Docker Hub | 認証と対象repositoryのpull/push権限をトークンのgrantで確認済み。イメージのアップロードは行っていない |

## 検証結果

- Python/R/Tosa全テスト: **252件成功、skip 0**。Rの全イベント、係数・contrast・交互作用、予測PSI、フィルタ/FDR、非推定可能な設計、収束失敗の扱いを含む。
- 合成BAMの全工程: Shiba Fisher / Shiba factorial / SnakeShiba Fisher / SnakeShiba factorial / scShiba / SnakeScShibaの**6通り成功**。bulkはStringTie、featureCounts、DESeq2、PCA、Excel、HTMLを含む。
- single-cell BAM/自己完結型（no_ref）CRAM: CB/UB、UMI重複除去、ライブラリをまたぐbarcode groupの集計、RI PSIが既知の期待値と一致。
- Excel: 混在/全空シート、値とヘッダー、freeze panes、出力失敗時の既存ファイル保護を確認。回帰モードのworkbookも生成確認。
- ブラウザ: 6通りの合成データ出力と公開データ出力で、外部ネットワークを遮断し、HTML表示、Plotlyズーム、SVG/PNG出力を確認。静的PNG/PDF/SVGの存在・ファイル形式も確認。
- `mkdocs build --strict`、workflowのactionlint、Python構文チェック成功。
- Docker統合CI: [`9da70b2`の全ジョブ成功](https://github.com/Sika-Zheng-Lab/Shiba/actions/runs/43503585120)。Shiba/MameShiba両候補のConda build・新規環境でのインストール/CLI/R起動テスト成功。公式channelの同一archiveと依存関係の修正metadataをローカルへミラーして検証した。

依存環境はPython 3.12.16、NumPy 2.2.6、pandas 2.3.3、SciPy 1.15.2、Tosa 1.0.0、Snakemake 9.13.7、R 4.5.3、glmmTMB 1.1.13、TMB 1.9.19、Matrix 1.7.5、DESeq2 1.50.2。CI artifactには実際のConda explicit spec、pip freeze、runtime.json、統合テスト出力を保存する（保持7日）。

## 公開データと測定範囲

| データ | 入力・範囲 | 実行結果 |
| --- | --- | --- |
| ERP133481 / SF3B1 | 既存TosaベンチマークのControl 2 + SF3B1 KD 2検体。GRCh38 chromosome 1:0–10 Mb、9,877,868 alignment records | Fisher全工程144.4秒 / peak PSS 933.5 MiB。回帰は集計済み入力から130.4秒 / 924.9 MiB |
| GSE246705 / Nanopore | 既存long-readベンチマークのcortex male 2 + female 2検体。GRCm38 chromosome 1:3–13 Mb、15,823 alignment records | Fisher全工程48.6秒 / 925.6 MiB。回帰は集計済み入力から36.3秒 / 919.8 MiB |
| 10x PBMC 1k v3 | 公開BAMのGRCh38 chromosome 1:0–10 Mb。filtered matrixに含まれる1,221 barcode、CB/UB付き764,545 alignment records | scShiba全工程22.5秒 / 198.3 MiB。barcodeを機械的に二分した実行検証であり、生物学的群間比較ではない |

10xデータは[1k PBMCs from a Healthy Donor (v3 chemistry)](https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0)、10x Genomics、2018、CC BY 4.0。取得元と抽出条件は`test/prepare_public_sc_data.py`に明記。bulkの抽出は`test/prepare_real_data.py`で再現できる。

短鎖データでは全8種類計2,245イベント（SE 393、FIVE 147、THREE 200、MXE 6、RI 216、MSE 176、AFE 796、ALE 311）を生成した。MXE 6件はすべてminus-strand。長鎖でも全8種類計230イベントを生成した。短鎖回帰はcoverage filter通過570件、長鎖は12件。実行成功と、すべてのイベントが推定可能であることは別であり、低coverage、untestable、方向不一致のstatusを保持する。

測定は`-p 2`、BLAS/OMP各1 thread。メモリはプロセスツリーの合計PSSを約1秒間隔で測定した。共有HPC上の単回測定であり、他の同時処理やキャッシュに影響される。上記は指定領域の検証であり、全ゲノム・大規模コホートの性能保証ではない。

## 回帰モデルの校正と制約

1,000イベント、2×2 design、各cell 8反復（32 samples）、depth 100、rho 0.05、seed 731の合成データを実バックエンドで検証した。処理80.7秒、監視プロセス込み86.4秒、peak PSS 862.0 MiB、untestable 0。帰無イベントの未補正p < 0.05の比率は8.0%（95%区間5.78–10.73%）、実現false discovery proportionは4.46%、powerは98.6%だった。

LRTは漸近近似であり、この条件でも未補正p値は名目5%より多く棄却された。単回の実現FDPを期待FDRの推定や普遍的なFDR制御の証明として扱わない。少数の生物学的反復、低coverage、収束しないモデルは特に注意が必要。長鎖の小領域・2反復/群は動作確認用であり、科学的校正の根拠にはしない。モデルはイベント間でdispersionを共有せず、single-cellへの独立反復を仮定した回帰には対応しない。

同条件の10,000イベント拡張でも成功した。処理744.7秒（監視込み761.7秒）、peak PSS 1085.3 MiB、untestable 0。未補正の帰無棄却率7.12%（95%区間6.42–7.87%）、実現FDP 3.80%、power 98.78%。この大きさでも名目水準を上回る未補正棄却が見られるため、漸近近似の制約として明記する。

## 検証記録

長期保存用の入力ID・抽出領域・checksum・依存バージョン・計測値・各statusを[`test/validation/release-2026-11-13.json`](test/validation/release-2026-11-13.json)に保存した。生ログ・出力は検証ホストの`/tmp/shiba-release-fixes/`にあり、一時保存である。Docker検証の出力は上記CIのartifactから取得できる。

## 再検証と配布手順

```bash
conda env create -f environment.yml
conda activate shiba
export PYTHONNOUSERSITE=1
python scripts/verify_runtime.py
python scripts/run_tests.py
python test/release_smoke.py /tmp/shiba-release-smoke
mkdocs build --strict
```

ブラウザ検証はPlaywright/Chromiumを追加導入して`python test/check_report_browser.py /tmp/shiba-release-smoke`を実行する。CIでは両Dockerfileで同じ手順を実行する。

Conda候補の生成・buildは、VERSIONを変更せず実施できる:

```bash
python scripts/prepare_conda_recipes.py --output /tmp/shiba-recipes
conda build /tmp/shiba-recipes/shiba --no-anaconda-upload --override-channels -c conda-forge -c bioconda
conda build /tmp/shiba-recipes/mameshiba --no-anaconda-upload --override-channels -c conda-forge -c bioconda
```

このローカル候補は現VERSIONを使うため、既存の正式版と混同して配布しない。

実際に公開する際は以下の順序にする:

1. 保留中のVERSIONとリリース日を確定し、`scripts/update_version.py`と`--check`でラベルを同期する。旧版に関する説明は維持する。
2. developの最終CI成功を確認し、develop → mainのPRを作成する。マージが本番公開のトリガー。
3. ワークフローがmainの対象SHA、VERSION、日付、既存タグを検証し、同じSHAのテスト・Docker全工程検証を実施する。
4. 検証済みイメージをversion tag、latestへ公開する。同じversionに別SHAのイメージがあれば上書きせず停止する。
5. タグを作成/再利用し、決定的ソースarchive・SHA256 manifest・image manifest・Shiba/MameShibaレシピをdraft Releaseへ登録してから公開する。ドキュメントは成功後に公開する。
6. 添付レシピをBiocondaの対象recipesへ提出し、Bioconda側のbuild/レビューと公開後の新規インストールを確認する。外部提出・本番公開は今回の準備作業では実施していない。

途中失敗時は同じSHAのworkflowを再実行する。既存の公開Releaseは成果物を維持する。mainが進んでいる場合やタグ/イメージが別SHAならpreflightが停止するため、古い候補でlatestを戻す操作はしない。
