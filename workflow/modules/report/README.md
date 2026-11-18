# report

One report of the run in `reports/`:

- `report.html`: one self-contained page (figures embedded)
- `report.docx`: Word, with a table of contents
- `report_notion.docx`: plain Word that pastes cleanly into Notion (headings
  1-3, short plain tables, PNG figures, captions as plain text)
- `report.md`: the Markdown they are made from

Built with `snake_core.report` (pandoc via `pypandoc_binary`; no Quarto).

## Sections

Summary; samples and their status; Cell Ranger metrics; preprocessing (cells
kept at each step, doublet, dsb and RNA QC figures per sample); TCR
clonotypes; the integrated object (metadata matches, metadata UMAPs and QC
plots); per population the cells per cluster, clustree, and for each shown
resolution the top markers, dot plots and marker plots; settings.

## Options (`report`)

| Option | Default | Meaning |
|---|---|---|
| `formats` | `[html, docx, notion]` | Formats to write |
| `title` | project name and version | Report title |
| `max_table_rows` | 40 | Rows shown per table (full tables stay in `results/`) |
| `top_markers` | 5 | Markers per cluster in the marker tables |
| `threshold_set` | `default` | Threshold set whose filtered markers are shown |
| `resolutions` | all | Resolutions shown per population, e.g. `{main: [0.07]}` |

## Rules

<!-- RULES:START -->
| Rule | Kind | What it does |
|---|---|---|
| `build_report` | rule | Write the run's report as Markdown and the formats in report.formats (HTML, DOCX, Notion DOCX). |
<!-- RULES:END -->
