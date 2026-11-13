"""Open generated reports offline and exercise Plotly zoom and SVG export."""
import argparse
import json
from pathlib import Path

from playwright.sync_api import sync_playwright


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    reports = sorted(args.directory.resolve().glob("*/plots/summary.html"))
    if not reports:
        raise RuntimeError("No reports found")
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True, args=["--no-sandbox"])
        results = []
        for report in reports:
            page = browser.new_page(accept_downloads=True)
            errors, network = [], []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.on("request", lambda request: network.append(request.url) if request.url.startswith(("http://", "https://")) else None)
            page.route("http://**/*", lambda route: route.abort())
            page.route("https://**/*", lambda route: route.abort())
            page.goto(report.as_uri())
            page.locator(".shiba-plot").first.scroll_into_view_if_needed()
            page.wait_for_function("document.querySelector('.shiba-plot').data !== undefined")
            result = page.evaluate("""async () => {
                const plot = document.querySelector('.shiba-plot');
                await Plotly.relayout(plot, {'xaxis.range': [0, 1]});
                const range = plot.layout.xaxis.range;
                const svg = await Plotly.toImage(plot, {format: 'svg'});
                const png = await Plotly.toImage(plot, {format: 'png'});
                return {range, svg: svg.startsWith('data:image/svg+xml'), png: png.startsWith('data:image/png')};
            }""")
            assert result["range"] == [0, 1] and result["svg"] and result["png"], result
            assert not errors, errors
            assert not network, network
            page.screenshot(path=str(report.parent / "browser-check.png"), full_page=False)
            for extension, signature in (("png", b"\x89PNG"), ("pdf", b"%PDF"), ("svg", b"<?xml")):
                exports = list(report.parent.rglob(f"*.{extension}"))
                assert exports, f"No {extension} exports for {report}"
                assert all(path.read_bytes().startswith(signature) for path in exports), extension
            results.append({"report": str(report), "offline": True, "zoom": True,
                            "svg_export": True, "png_export": True, "static_formats": ["png", "pdf", "svg"]})
            page.close()
        browser.close()
    (args.directory / "browser-validation.json").write_text(json.dumps(results, indent=2))
    print(json.dumps(results))


if __name__ == "__main__":
    main()
