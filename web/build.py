#!/usr/bin/env python3
"""Assemble the web app.

  python3 build.py [--dev]

writes  index.html          -- development page that loads vendor/ and src/ files separately
        dist/index.html     -- the single self-contained file (open it directly, host it anywhere)
        dist/artifact.html  -- the same page without the document wrapper, for publishing as a claude.ai artifact

--dev   also offers the examples that are not elliptic curves (genus 2, a nodal cubic, a cubic not in
        Weierstrass form), which exercise the error paths; the pages then define window.EC3D_DEV.
"""
import argparse, os
HERE = os.path.dirname(os.path.abspath(__file__))
rd = lambda p: open(os.path.join(HERE, p), encoding="utf-8").read()

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dev", action="store_true", help="include the examples that are not elliptic curves")
args = ap.parse_args()

import base64, re
template, body, css = rd("src/index.template.html"), rd("src/body.html"), rd("src/app.css")
scripts = ["vendor/three.min.js", "vendor/OrbitControls.js", "vendor/katex/katex.min.js", "src/colormaps.js", "src/curves-data.js", "src/ec3d-math.js", "src/app.js"]

def katex_css_inline():
    """katex.min.css with the vendored woff2 fonts embedded as data URIs (the woff/ttf fallbacks dropped)."""
    k = rd("vendor/katex/katex.min.css")
    def font(m):
        name = m.group(1)
        path = os.path.join(HERE, "vendor/katex/fonts", name + ".woff2")
        if not os.path.exists(path):
            return m.group(0)
        data = base64.b64encode(open(path, "rb").read()).decode()
        return f"url(data:font/woff2;base64,{data}) format(\"woff2\")"
    k = re.sub(r'url\(fonts/(KaTeX_[\w-]+)\.woff2\) format\("woff2"\)', font, k)
    return re.sub(r',url\(fonts/KaTeX_[\w-]+\.(?:woff|ttf)\) format\("(?:woff|truetype)"\)', '', k)
for s in scripts:
    assert "</script" not in rd(s).lower(), s + " contains </script"
flags = "<script>window.EC3D_DEV = true;</script>\n" if args.dev else ""
dev = template.replace("<!--STYLES-->", '<link rel="stylesheet" href="vendor/katex/katex.min.css">\n<link rel="stylesheet" href="src/app.css">').replace("<!--BODY-->", body) \
              .replace("<!--SCRIPTS-->", flags + "\n".join(f'<script src="{s}"></script>' for s in scripts))
inline_scripts = flags + "\n".join(f"<script>\n{rd(s)}\n</script>" for s in scripts)
styles = f"<style>\n{katex_css_inline()}\n</style>\n<style>\n{css}\n</style>"
dist = template.replace("<!--STYLES-->", styles).replace("<!--BODY-->", body).replace("<!--SCRIPTS-->", inline_scripts)
artifact = f"<title>Elliptic Curve 3D</title>\n{styles}\n{body}\n{inline_scripts}\n"
os.makedirs(os.path.join(HERE, "dist"), exist_ok=True)
for name, content in [("index.html", dev), ("dist/index.html", dist), ("dist/artifact.html", artifact)]:
    open(os.path.join(HERE, name), "w", encoding="utf-8").write(content)
    print(f"wrote {name}: {len(content.encode()):,} bytes" + (" (dev examples included)" if args.dev else ""))
