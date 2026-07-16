import markdown, re
from weasyprint import HTML

md = open("BIOLOGY_BACKGROUND.md", encoding="utf-8").read()
html_body = markdown.markdown(md, extensions=["fenced_code", "tables", "codehilite", "toc"])
# render the hourglass emoji reliably as a red badge
html_body = html_body.replace("⏳", '<span class="flag">&#9873; &gt;4h</span>').replace(" &gt;4h &gt;4h", " &gt;4h")

css = """
@page { size: A4; margin: 1.6cm 1.6cm 1.8cm 1.6cm;
        @bottom-center { content: counter(page) " / " counter(pages);
                         font-size: 8pt; color: #888; } }
* { box-sizing: border-box; }
body { font-family: -apple-system, 'Helvetica Neue', Arial, sans-serif;
       font-size: 10pt; line-height: 1.5; color: #1a1a1a; }
h1 { font-size: 20pt; color: #0b3d5c; border-bottom: 3px solid #0b3d5c;
     padding-bottom: 6px; margin-top: 0; }
h2 { font-size: 14pt; color: #0b3d5c; border-bottom: 1px solid #ccc;
     padding-bottom: 3px; margin-top: 22px; page-break-after: avoid; }
h3 { font-size: 11.5pt; color: #14567a; margin-top: 18px; page-break-after: avoid; }
p, li { orphans: 2; widows: 2; }
a { color: #14567a; text-decoration: none; word-break: break-all; }
code { font-family: 'SFMono-Regular', Menlo, Consolas, monospace; font-size: 8.6pt;
       background: #f2f4f6; padding: 1px 4px; border-radius: 3px; color: #b8003a; }
pre { background: #1e2733; color: #e6edf3; padding: 10px 12px; border-radius: 6px;
      font-size: 8.2pt; line-height: 1.4; overflow-wrap: break-word;
      white-space: pre-wrap; page-break-inside: avoid; }
pre code { background: none; color: #e6edf3; padding: 0; }
table { border-collapse: collapse; width: 100%; font-size: 8.8pt; margin: 10px 0;
        page-break-inside: avoid; }
th { background: #0b3d5c; color: #fff; text-align: left; padding: 5px 7px; }
td { border: 1px solid #d0d7de; padding: 5px 7px; vertical-align: top; }
tr:nth-child(even) td { background: #f6f8fa; }
hr { border: none; border-top: 1px solid #d0d7de; margin: 20px 0; }
.flag { color: #c0392b; font-weight: 700; white-space: nowrap; }
blockquote { border-left: 3px solid #14567a; margin-left: 0; padding-left: 12px; color: #444; }
"""

html = f"<html><head><meta charset='utf-8'><style>{css}</style></head><body>{html_body}</body></html>"
HTML(string=html).write_pdf("BIOLOGY_BACKGROUND.pdf")
print("wrote BIOLOGY_BACKGROUND.pdf")
