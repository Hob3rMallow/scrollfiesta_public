#!/usr/bin/env python3
"""make_gallery_pdf.py -- tile per-mesh PNGs into a labelled multi-page PDF
contact sheet for spotting garbage sheets at a glance.

PNGs are sorted by filename (= cube id z<.>_y<.>_x<.>, i.e. spatial z,y,x order)
so neighbours sit together and a bad sheet stands out against its neighbours.
Each tile is labelled with its cube id. PIL only (no reportlab/img2pdf).

Usage: python make_gallery_pdf.py <png_dir> <out.pdf> [--cols N] [--rows N] [--tile PX]
"""
import sys, os, glob
from PIL import Image, ImageDraw, ImageFont

def main():
    if len(sys.argv) < 3:
        print("usage: make_gallery_pdf.py <png_dir> <out.pdf> [--cols N] [--rows N] [--tile PX]")
        return 1
    png_dir, out_pdf = sys.argv[1], sys.argv[2]
    cols, rows, tile = 6, 8, 360
    a = sys.argv[3:]
    for i, t in enumerate(a):
        if t == "--cols": cols = int(a[i+1])
        elif t == "--rows": rows = int(a[i+1])
        elif t == "--tile": tile = int(a[i+1])

    pngs = sorted(glob.glob(os.path.join(png_dir, "*.png")))
    if not pngs:
        print(f"no PNGs in {png_dir}"); return 1
    print(f"{len(pngs)} tiles, {cols}x{rows}/page, tile={tile}px")

    try:
        font = ImageFont.truetype("arial.ttf", 16)
    except Exception:
        font = ImageFont.load_default()

    pad, lbl_h = 6, 20
    cell_w = tile + pad
    cell_h = tile + lbl_h + pad
    page_w = cols * cell_w + pad
    page_h = rows * cell_h + pad
    per_page = cols * rows

    pages = []
    for p0 in range(0, len(pngs), per_page):
        chunk = pngs[p0:p0 + per_page]
        page = Image.new("RGB", (page_w, page_h), (255, 255, 255))
        draw = ImageDraw.Draw(page)
        for k, path in enumerate(chunk):
            r, c = divmod(k, cols)
            x = pad + c * cell_w
            y = pad + r * cell_h
            try:
                im = Image.open(path).convert("RGB")
                im.thumbnail((tile, tile), Image.LANCZOS)
                page.paste(im, (x + (tile - im.width)//2, y + (tile - im.height)//2))
            except Exception as e:
                draw.rectangle([x, y, x+tile, y+tile], outline=(200,0,0))
            cube = os.path.splitext(os.path.basename(path))[0]
            draw.text((x + 2, y + tile + 2), cube, fill=(0, 0, 0), font=font)
        pages.append(page)
        print(f"  page {len(pages)} ({len(chunk)} tiles)")

    pages[0].save(out_pdf, save_all=True, append_images=pages[1:])
    sz = os.path.getsize(out_pdf) / (1024*1024)
    print(f"wrote {out_pdf}  ({len(pages)} pages, {sz:.1f} MB)")
    return 0

if __name__ == "__main__":
    sys.exit(main())
