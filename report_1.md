# Отчет по первому филогенетическому практикуму
## Создание филогенетических деревьев в R и Python

### Начнем с Python
1. Изначально прочитаем дерево с помощью модуля *requests**. После приступаем к отрисовке дерева с псевдографикой *draw_ascii* в консоли. Далее отрисуем дерево и сохраним картинку:
```python
# чтение дерева
raw_tree = StringIO(requests.get('https://www.jasondavies.com/tree-of-life/life.txt').text)
tree1 = Phylo.read(raw_tree, "newick")

# отрисовка дерева в консоли
Phylo.draw_ascii(tree1)
```

### R
