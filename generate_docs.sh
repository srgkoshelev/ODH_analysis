echo 'Generating docs'
pdoc --html ODH_analysis -o docs
mv docs/ODH_analysis/* docs
rmdir docs/ODH_analysis
echo 'Docs generated'
git add -A
git commit -m 'Generated docs.'
git push
echo 'Committed and pushed upstream.'
