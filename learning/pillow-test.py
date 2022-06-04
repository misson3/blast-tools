# Jun04, 2022, ms
# pillow-test.py

from PIL import Image, ImageDraw

img = Image.new("RGB", (300, 100), "White")
draw = ImageDraw.Draw(img)

draw.line([(10, 90), (290, 10)], fill="Blue")

img.show()
