# 05_generate_report.R
# ---------------------
# Render RMarkdown report

rmarkdown::render("reports/analysis_report.Rmd", output_format = "html_document")
