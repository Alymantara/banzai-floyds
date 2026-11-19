"""Build quickstart_reduction.ipynb from quickstart_reduction.py (cells separated by '# %%', markdown cells by
'# %% [markdown]' with every line commented), so the notebook and the script that is tested stay the same.

usage: python notebooks/make_notebook.py"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def cells(source: str) -> list:
    out = []
    for block in source.split('# %%')[1:]:
        markdown = block.startswith(' [markdown]')
        lines = block.split('\n')[1:]
        while lines and not lines[-1].strip():
            lines.pop()
        if markdown:
            lines = [line[2:] if line.startswith('# ') else line.lstrip('#') for line in lines]
            out.append({'cell_type': 'markdown', 'metadata': {}, 'source': [line + '\n' for line in lines[:-1]] +
                        lines[-1:]})
        else:
            out.append({'cell_type': 'code', 'metadata': {}, 'execution_count': None, 'outputs': [],
                        'source': [line + '\n' for line in lines[:-1]] + lines[-1:]})
    return out


if __name__ == '__main__':
    source = open(os.path.join(HERE, 'quickstart_reduction.py')).read()
    notebook = {'cells': cells(source), 'nbformat': 4, 'nbformat_minor': 5,
                'metadata': {'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
                             'language_info': {'name': 'python'}}}
    json.dump(notebook, open(os.path.join(HERE, 'quickstart_reduction.ipynb'), 'w'), indent=1)
    print(len(notebook['cells']), 'cells written to quickstart_reduction.ipynb')
