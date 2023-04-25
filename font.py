from PIL import Image, ImageDraw, ImageFont

# Define the font names and sample text
font_names = ["Arial", "Balto", "Courier New", "Droid Sans", "Droid Serif",
         "Droid Sans Mono", "Gravitas One", "Old Standard TT", "Open Sans",
         "Overpass", "PT Sans Narrow", "Raleway", "Times New Roman"]
text = "Hello, world!"

# Set up the image and draw object
width, height = 800, 600
image = Image.new("RGB", (width, height), color="white")
draw = ImageDraw.Draw(image)

# Loop over the font names and draw the text in each font
y = 50
for font_name in font_names:
    font = ImageFont.truetype(font_name + ".ttf", size=36)
    draw.text((50, y), text, font=font, fill="black")
    y += 50

# Show the resulting image
image.show()