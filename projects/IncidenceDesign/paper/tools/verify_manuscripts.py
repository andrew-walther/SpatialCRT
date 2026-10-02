#!/usr/bin/env python3
"""Check current manuscript tables against verified aggregate outputs.

Author: Codex (reviewed by Andrew Walther). Created: 2026-10-02.
Dependencies: Python standard library.

Run from any directory: python3 paper/tools/verify_manuscripts.py
Only aggregate CSVs and document sources are read. This does not assert that
all prose claims are correct: independent scientific review covers their scope.
"""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CHAPTER = ROOT / 'paper/dissertation_chapter/Dissertation_Chapter.qmd'
MAIN = ROOT / 'paper/ctj_manuscript/CTJ_Manuscript.tex'
MAIN_EXHIBITS = MAIN.with_name('CTJ_Exhibits.tex')
READING = MAIN.with_name('CTJ_Reading.tex')
SUPPLEMENT = ROOT / 'paper/ctj_manuscript/Supplementary_Information.tex'
EXHIBITS = ROOT / 'application/results/real_sud_rev_20261002/exhibits'
LABELS = ['Graph Checkerboard', 'High Incidence Focus', 'Plain saturation',
          'Isolation Buffer', 'Spatial blocking', 'Balanced Quartiles',
          'Balanced Halves', 'Guided saturation', 'SRS']


def csv_rows(path):
    """Read authorized aggregate records as strings, preserving printed precision."""
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def table(source, label):
    """Extract the tabular body belonging to an explicit LaTeX table label."""
    position = source.index('\\label{' + label + '}')
    start = source.index('\\begin{tabular}', position)
    start = source.index('\n', start) + 1
    end = source.index('\\end{tabular}', start)
    return source[start:end]


def numerical_rows(body):
    """Return decimal cells by row label; headers without numbers are excluded."""
    result = []
    for line in body.splitlines():
        if '&' not in line or line.lstrip().startswith('%'):
            continue
        cells = line.split('&')
        if cells[0].strip() in ('', 'Design', 'Sensitivity') or cells[0].lstrip().startswith('\\'):
            continue
        numbers = re.findall(r'(?<![\w.])-?\d+(?:\.\d+)?', '&'.join(cells[1:]))
        if numbers:
            result.append((cells[0].strip(), numbers))
    return result


def require_equal(actual, expected, context):
    """Reject a wrong printed value or order with its source-table context."""
    if actual != expected:
        raise AssertionError(f'{context}: {actual!r} != {expected!r}')


def verify():
    """Check application/grid tables, chapter derivation, assets and citations."""
    sources = {p: p.read_text() for p in (CHAPTER, MAIN, SUPPLEMENT)}
    # The master owns prose; both layouts use these exact exhibit definitions.
    master = sources[MAIN]
    shared_exhibits = MAIN_EXHIBITS.read_text()
    reading = READING.read_text()
    legends = MAIN.with_name('Figure_Legends.md').read_text()
    figure_blocks = re.findall(r'\\begin\{figure\*\}.*?\\end\{figure\*\}',
                               shared_exhibits, re.S)
    require_equal(len(figure_blocks), 4, 'Four main figures')
    for index, block in enumerate(figure_blocks, 1):
        caption = re.search(r'\\caption\{(.*?)\}\\label\{fig:', block, re.S).group(1)
        caption = caption.replace('Table~\\ref{tab:grid-srs}', 'Table 2')
        if f'Figure {index}. {caption}' not in legends:
            raise AssertionError(f'Figure {index}: separate legend differs from master')
    require_equal(re.search(r'\\newcommand\{\\CTJTitle\}\{([^}]+)\}', master).group(1),
                  'Sampling Design for Spatial Cluster Randomized Trials Under Heterogeneous Incidence',
                  'Author-approved restored title')
    author_block = master.split('\\newcommand{\\CTJAuthorLine}', 1)[1].split('\n', 1)[0]
    authors = ['Andrew Walther', 'Ashkan Habib', 'Ross Joseph Simpson, Jr.', 'Feng-Chang Lin']
    positions = [author_block.index(author) for author in authors]
    require_equal(positions, sorted(positions), 'Author-approved author order')
    require_equal(reading.count('\\input{CTJ_Manuscript.tex}'), 1,
                  'Reading preview uses the master manuscript')
    require_equal(reading.count('\\def\\CTJReadingVersion{1}'), 1,
                  'Reading preview selects the alternate layout')
    if re.search(r'\\(?:section|caption|includegraphics)\b', reading):
        raise AssertionError('Reading wrapper duplicates manuscript content')
    require_equal(master.count('\\input{CTJ_Exhibits.tex}'), 1,
                  'Master loads the shared exhibit definitions')
    for name in ('CTJDesignTable', 'CTJGridTable', 'CTJMapFigure',
                 'CTJGridFigure', 'CTJBiasFigure', 'CTJMeanFigure'):
        calls = [m.start() for m in re.finditer(r'\\' + name + r'\b', master)]
        require_equal(len(calls), 2, name + ': one call in each layout branch')
        bibliography = master.index('\\bibliographystyle')
        if not calls[0] < bibliography < calls[1]:
            raise AssertionError(name + ': embedded/end-positioned placement')
        require_equal(shared_exhibits.count('\\newcommand{\\' + name + '}'), 1,
                      name + ': one shared definition')
    sources[MAIN] += '\n' + shared_exhibits
    chapter, main, supplement = (sources[p] for p in (CHAPTER, MAIN, SUPPLEMENT))
    for path, source in ((CHAPTER, chapter), (SUPPLEMENT, supplement)):
        author_line = next(line for line in source.splitlines() if '\\author{' in line)
        positions = [author_line.index(author) for author in authors]
        require_equal(positions, sorted(positions), path.name + ': author order')
    annual = csv_rows(EXHIBITS / 'yearly_primary_design_means.csv')
    yearly = {(r['Year'], r['Regime'], int(r['Design_ID'])): r for r in annual}
    grid = csv_rows(ROOT / 'results/srs_benchmark/manuscript_regime_means.csv')
    grid = {(r['Design'], r['Spillover_Type']): r for r in grid}
    expected_grid = []
    for design_id in range(1, 10):
        cells = []
        for regime in ('both', 'control_only'):
            row = grid[(f'Design {design_id}', regime)]
            cells.extend([f"{float(row['MSE']):.4f}", f"{float(row['Coverage']):.3f}"])
        expected_grid.append(cells)
    expected_year = []
    for year in ('2018', '2019', '2020', '2021'):
        cells = [f"{float(yearly[(year, regime, design)]['Mean_MSE']):.5f}"
                 for regime, design in [('both', 9), ('both', 8),
                                        ('control_only', 9), ('control_only', 3),
                                        ('control_only', 8)]]
        expected_year.append(cells)
    for path, source in sources.items():
        require_equal([n for _, n in numerical_rows(table(source, 'tab:grid-srs'))],
                      expected_grid, f'{path.name}: grid benchmark')
        if path != MAIN:
            require_equal([n for _, n in numerical_rows(table(source, 'tab:nc-year'))],
                          expected_year, f'{path.name}: primary annual MSE')
    # The main now reports annual performance graphically; validate all 72 ratios
    # against the same primary aggregates used by the chapter/SI tables.
    compact = ROOT / 'results/manuscript_exhibit_revision_20261002'
    ratios = csv_rows(compact / 'annual_srs_ratios.csv')
    require_equal(len(ratios), 72, 'Compact annual figure cell count')
    for row in ratios:
        key = (row['Year'], row['Regime'], int(row['Design_ID']))
        expected = float(yearly[key]['Mean_MSE']) / float(yearly[(key[0], key[1], 9)]['Mean_MSE'])
        if abs(float(row['Ratio']) - expected) > 1e-10:
            raise AssertionError(f'Compact annual figure: wrong ratio {key}')
    require_equal(len(re.findall(r' & ', table(main, 'tab:allocation-rules'))), 10,
                  'Allocation-rule table: header plus nine strategies')
    for year in ('2018', '2019', '2020', '2021'):
        expected = []
        for regime in ('both', 'control_only'):
            for design in range(1, 10):
                row = yearly[(year, regime, design)]
                expected.append([f"{float(row[column])*scale:.{digits}f}"
                                 for column, digits, scale in
                                 [('Mean_MSE', 5, 1), ('Bias', 4, 1), ('Coverage', 3, 1),
                                  ('Mean_Treated', 2, 1), ('Mean_Population_Share', 1, 100)]])
        actual = numerical_rows(table(chapter, 'tab:nc-full-' + year))
        require_equal([n for _, n in actual], expected, f'{year}: full annual performance')
        require_equal([name for name, _ in actual], LABELS * 2, f'{year}: annual design order')
    budget_expected = []
    for design in range(1, 10):
        cells = []
        for year in ('2018', '2019', '2020', '2021'):
            row = yearly[(year, 'both', design)]
            cells.extend([f"{float(row['Mean_Treated']):.2f}",
                          f"{100*float(row['Mean_Population_Share']):.1f}"])
        budget_expected.append(cells)
    require_equal([n for _, n in numerical_rows(table(chapter, 'tab:nc-budgets'))],
                  budget_expected, 'Primary annual cluster/population budgets')
    sensitivity = csv_rows(EXHIBITS / 'yearly_sensitivity_means.csv')
    matched = {(r['Model'], r['Neighbor'], r['Summary'], r['Year'], r['Regime']): r
               for r in sensitivity if r['Design_ID'] == '8'}
    expected_sensitivity = []
    for model, neighbor, summary in [('baseline_sensitivity', 'queen', 'mean_rank'),
                                     ('education', 'rook', 'mean_rank'),
                                     ('education', 'queen', 'mean_rate'),
                                     ('education', 'queen', 'pooled_rate')]:
        for regime in ('both', 'control_only'):
            expected_sensitivity.append([
                f"{float(matched[(model, neighbor, summary, year, regime)]['Ratio_of_Matched_Mean_MSE']):.3f}"
                for year in ('2018', '2019', '2020', '2021')])
    require_equal([n for _, n in numerical_rows(table(chapter, 'tab:nc-summary-sens'))],
                  expected_sensitivity, '32 parameter-matched sensitivity ratio cells')
    # Check every supplement table body against its long-form chapter source.
    table_labels = re.findall(r'\\label\{(tab:[^}]+)\}', chapter)
    for label in table_labels:
        def normalized(body):
            """Ignore comments/whitespace, retaining every cell and table command."""
            body = body.replace('Chapter 2', 'the earlier small-grid study')
            return re.sub(r'\s+', '', re.sub(r'(?<!\\)%[^\n]*', '', body))
        require_equal(normalized(table(chapter, label)), normalized(table(supplement, label)),
                      f'Supplement derivation: {label}')
    bib = (ROOT / 'paper/SpatialCRT_IncidenceDesign.bib').read_text()
    keys = set(re.findall(r'@\w+\s*\{([^,]+),', bib))
    for path, source in sources.items():
        assets = re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', source)
        for asset in assets:
            if not (path.parent / asset).is_file():
                raise AssertionError(f'{path.name}: missing figure {asset}')
        reference_source = re.sub(r'(?<!\\)%[^\n]*', '', source)
        reference_source = re.sub(r'`[^`]*`', '', reference_source)
        refs = re.findall(r'\\ref\{([^}]+)\}', reference_source)
        labels = set(re.findall(r'\\label\{([^}]+)\}', source))
        labels.update(re.findall(r'\{#([^}]+)\}', source))
        labels.update(re.findall(r'\\applabel\{[^}]+\}\{([^}]+)\}', source))
        require_equal(set(refs) - labels, set(), f'{path.name}: unresolved source references')
        cites = []
        for block in re.findall(r'\\cite\w*(?:\[[^]]*\])*\{([^}]+)\}', source):
            cites.extend(key.strip() for key in block.split(','))
        if path == CHAPTER:
            prose = source.split('---', 2)[-1]
            cites.extend(re.findall(r'(?<![\w])@([\w:.-]+)', prose))
        require_equal(set(cites) - keys, set(), f'{path.name}: unknown bibliography keys')
    require_equal(len(re.findall(r'\\begin\{(?:table|figure)\*?\}', main)), 6, 'CTJ exhibit cap')
    for source in (main, supplement):
        if re.search(r'Chapter 2|Project 1|NEEDS-AUTHOR-CONFIRMATION|\bDIM\b', source):
            raise AssertionError('Stale/non-standalone CTJ wording')
    # A deliberate erroneous value must fail the same numerical comparison.
    wrong = [list(row) for row in expected_grid]
    wrong[0][0] = '99.9999'
    try:
        require_equal(wrong, expected_grid, 'Wrong-value fixture')
    except AssertionError:
        pass
    else:
        raise AssertionError('Wrong-value fixture was silently accepted')
    print('PASS: 36 grid cells in all three sources; 20 annual MSE cells in chapter/SI;')
    print('72 matched annual figure ratios and all-nine-design allocation-rule table;')
    print('360 annual performance cells, 72 budget cells, 32 sensitivity cells, all ' + str(len(table_labels)) +
          ' supplement table bodies; figures, references, citations and six-exhibit cap.')
    print('PASS: deliberate incorrect-number fixture rejected. Scientific prose reviewed separately.')
    print('PASS: reading/submission share one master and six exhibit definitions; '
          'embedded/end-positioned calls checked.')


if __name__ == '__main__':
    verify()
