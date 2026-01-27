# install.packages("magick")  # if not already installed
library(magick)

# read the original image
img <- image_read("Fig 7.tif")

# add 0.25-inch white border on all sides
# (adjust values: "25x25" = 25 pixels; "50x50" = thicker)
bordered <- image_border(img, color = "white", geometry = "50x50")

# save it again (lossless)
image_write(bordered, path = "Fig 7.tiff", format = "tiff")
