// Script: test_companion.mjs
// Purpose: Exercise the saved companion's JavaScript with a minimal DOM fixture.
// Author: Andrew Walther
// Created: 2026-10-02
// Dependencies: Node built-ins only. This is not a visual browser review.
import fs from 'node:fs';
import vm from 'node:vm';
import assert from 'node:assert/strict';

const htmlPath = new URL('../report/real_sud_companion.html', import.meta.url);
const html = fs.readFileSync(htmlPath, 'utf8');
const scripts = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(m => m[1]);
assert.equal(scripts.length, 2);
const elements = {};

// Provide only the DOM operations used by the actual saved page. Unexpected
// operations fail; callbacks execute unchanged rather than being reimplemented.
class Element {
  constructor(tag) { this.tag = tag; this.children = []; this.value = ''; this.attrs = {}; }
  set id(value) { this._id = value; elements[value] = this; }
  get id() { return this._id; }
  append(child) {
    this.children.push(child);
    if (this.tag === 'select' && this.children.length === 1) this.value = child.value;
  }
  replaceChildren(...children) { this.children = children; }
  setAttribute(key, value) { this.attrs[key] = String(value); }
}
for (const id of ['filters', 'metric', 'dataset', 'rows', 'chart', 'resultStatus']) {
  const e = new Element(['metric', 'dataset'].includes(id) ? 'select' : 'div'); e.id = id;
}
elements.metric.value = 'Mean_MSE'; elements.dataset.value = 'main';
const context = vm.createContext({document: {
  getElementById: id => elements[id],
  createElement: tag => new Element(tag),
  createElementNS: (_ns, tag) => new Element(tag)
}});
for (const script of scripts) new vm.Script(script).runInContext(context);
assert.equal(elements.rows.children.length, 9, 'Initial primary setting must display nine designs');
assert.equal(elements.filters.children.length, 7);
const bars = elements.chart.children.filter(e => e.tag === 'rect');
assert.equal(bars.length, 9);
assert.ok(Number(bars[0].attrs.width) !== Number(bars[8].attrs.width), 'Bars must reflect differing MSEs');
assert.ok(bars.every(e => Number.isFinite(Number(e.attrs.width))));

// Check a planned summary-only cell, an unplanned combination, and refined data.
elements.Summary.value = 'mean_rate'; elements.Summary.onchange();
assert.equal(elements.rows.children.length, 1);
elements.Summary.value = 'mean_rank'; elements.Model.value = 'baseline_sensitivity';
elements.Neighbor.value = 'rook'; elements.Neighbor.onchange();
assert.equal(elements.rows.children.length, 0, 'Unplanned crossed sensitivities cannot fabricate rows');
elements.dataset.value = 'tail'; elements.dataset.onchange();
assert.equal(elements.rows.children.length, 9);
assert.equal(elements.Gamma.value, '0.8');
elements.Year.value = '2021'; elements.Year.onchange();
assert.equal(elements.rows.children.length, 9);
elements.metric.value = 'Q90_Estimated'; elements.metric.onchange();
assert.equal(elements.chart.children.filter(e => e.tag === 'rect').length, 9);
assert.equal(elements.chart.children.filter(e => e.tag === 'line').length, 0,
  'Mean-MSE MC intervals must not be presented as q90 uncertainty');
for (const match of html.matchAll(/<img src="([^"]+)"/g)) {
  assert.ok(fs.existsSync(new URL(match[1], htmlPath)), `Missing map/figure: ${match[1]}`);
}
console.log('PASS: actual companion scripts load, filter main/refined evidence, respect planned cells, scale values and avoid false tail intervals.');
console.log('This DOM fixture verifies behavior; visual browser rendering remains unverified.');
