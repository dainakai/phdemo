# Delete directories inside 20250121 except 00_img
find ./20250121 -mindepth 1 -maxdepth 1 -type d ! -name '00_img' -exec rm -rf {} +

# adjustted_image.png, after_BA.png, before_BA.png があれば削除
rm -f ./adjusted_image.png
rm -f ./after_BA.png
rm -f ./before_BA.png