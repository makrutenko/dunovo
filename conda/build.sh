
# Download submodules and move them to lib.
get_submodule () {
  read name version owner repo hash <<< "$@"
  wget --no-check-certificate "https://github.com/$owner/$repo/archive/v$version.tar.gz"
  downloaded_hash=`sha256sum "v$version.tar.gz" | tr -s ' ' | cut -d ' ' -f 1`
  if ! [ "$hash" == "$downloaded_hash" ]; then
    echo "Error: Hash does not match!" >&2
    return 1
  fi
  tar -zxvpf "v$version.tar.gz"
  rm -rf "$name"
  if [ "$name" == kalign ]; then
    mv "$repo-$version" "$name"
  else
    mv "$repo-$version" "$PREFIX/lib/$name"
  fi
  rm "v$version.tar.gz"
}
get_submodule kalign  0.3.0       makovalab-psu kalign-dunovo c0ef2de4a958aed47311ea86591debb3bada871143a6923bc6338ab2f99f2d5b
get_submodule utillib 0.1.1-alpha NickSto       utillib       961f5a3481d1c0dbe00c258b9df2de541e5605f3de4ff25bcc3cec22922e7c06
get_submodule ET      0.3         NickSto       ET            6b757b3ab3634b949f78692ac0db72e9264a1dceaf7449e921091d3ad0012eea
get_submodule bfx     0.4.0       NickSto       bfx           684f1f7bc9a8767bb1703addad1b669246723ff08915d404cecbbc1af4d3b3b3

# Compile binaries and move them to lib.
make
mv *.so "$PREFIX/lib"
mv kalign "$PREFIX/lib"

# Move scripts to lib and link to them from bin.
mkdir "$PREFIX/bin"
for script in *.awk *.sh *.py; do
  mv "$script" "$PREFIX/lib"
  ln -s "../lib/$script" "$PREFIX/bin/$script"
done
# Handle special cases.
mv utils/precheck.py "$PREFIX/lib"
ln -s ../lib/precheck.py "$PREFIX/bin"
ln -s ../lib/bfx/trimmer.py "$PREFIX/bin"
mv VERSION "$PREFIX/lib"
