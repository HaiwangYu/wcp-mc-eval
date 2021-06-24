# wcp-mc-eval

```
find . -type f -name "*.dot" | xargs dot -Tpng -O
rename -e 's/.dot.png/.png/' *.dot.png
```
