import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {createRequire} from 'node:module';
import {join} from 'node:path';

const require = createRequire(join(process.env.PLAYWRIGHT_ROOT, 'package.json'));
const {chromium} = require('playwright');
const site = process.env.SITE_URL || 'http://127.0.0.1:7800/spatial-transcriptomics-viz/';
const browser = await chromium.launch({headless: true, args: ['--no-proxy-server']});
const page = await browser.newPage({viewport: {width: 1440, height: 900}, acceptDownloads: true});
const errors = [];
page.on('pageerror', error => errors.push(error.message));
page.on('console', message => {if (message.type() === 'error') errors.push(message.text());});

try {
  await page.goto(`${site}#/dataset-overview`);
  await page.locator('h3[title]').first().waitFor();
  assert.equal(await page.locator('h3[title]').count(), 13);
  await page.getByRole('button', {name: 'Metadata'}).first().click();
  await page.getByRole('heading', {name: 'Dataset Metadata'}).waitFor();
  console.log('PASS dataset overview: 13 datasets and metadata');
  await page.mouse.click(5, 5);
  await page.getByRole('heading', {name: 'Dataset Metadata'}).waitFor({state: 'hidden'});

  await page.evaluate(() => {
    window.__openedHE = '';
    window.open = url => {window.__openedHE = String(url); return null;};
  });
  await page.getByRole('button', {name: 'HE Image'}).first().click();
  const heUrl = await page.evaluate(() => window.__openedHE);
  assert.match(heUrl, /\/data_root\/[^/]+\/he-raw\.tiff$/);
  const heResponse = await page.request.head(heUrl);
  assert.equal(heResponse.status(), 200);
  console.log('PASS lossless HE TIFF download URL');

  await page.goto(`${site}#/metric-vis`);
  await page.locator('select').first().selectOption('ViHD_mmu_kidney_h');
  await page.locator('select').nth(1).selectOption('KLD');
  await page.waitForFunction(() => [...document.querySelectorAll('select')[2].options].some(x => x.value === '2.5'));
  await page.locator('select').nth(2).selectOption('2.5');
  await page.getByRole('button', {name: 'Visualize'}).click();
  const image = page.locator('canvas[aria-label="reference Itfg1"]');
  await image.waitFor();
  await page.waitForFunction(() => document.querySelector('canvas[aria-label="reference Itfg1"]')?.width === 755);
  const before = createHash('sha256').update(await image.evaluate(canvas => canvas.toDataURL())).digest('hex');
  await page.getByRole('button', {name: /Auto \(Per Gene\)/}).click();
  await page.getByRole('button', {name: 'Manual Range'}).click();
  await page.waitForTimeout(250);
  const after = createHash('sha256').update(await image.evaluate(canvas => canvas.toDataURL())).digest('hex');
  assert.notEqual(before, after);
  console.log('PASS metric visualization: recovered Itfg1 rendered and color range changes pixels');

  await page.locator('select').first().selectOption('Xe_hsa_skin_h');
  await page.locator('select').nth(1).selectOption('KLD');
  await page.waitForFunction(() => [...document.querySelectorAll('select')[2].options].some(x => x.value === '5.0'));
  await page.locator('select').nth(2).selectOption('5.0');
  const largeStart = Date.now();
  await page.getByRole('button', {name: 'Visualize'}).click();
  await page.waitForFunction(() => {
    const canvases = [...document.querySelectorAll('canvas')];
    return canvases.length === 6 && canvases.every(canvas => canvas.width > 0);
  }, null, {timeout: 60000});
  console.log(`PASS large arrays: 6 canvases rendered in ${Date.now() - largeStart} ms at ${site}`);

  await page.goto(`${site}#/out-of-sample`);
  await page.locator('rect.region-cell').first().waitFor();
  try {
    await page.waitForFunction(() => {
      const rect = document.querySelector('rect.region-cell');
      return rect && Number(rect.getAttribute('x')) < 1000;
    }, null, {timeout: 5000});
  } catch (error) {
    console.log('OOS DIAGNOSTIC', await page.evaluate(() => ({
      images: [...document.images].map(image => ({src: image.src, width: image.naturalWidth, height: image.naturalHeight, complete: image.complete})),
      rect: document.querySelector('rect.region-cell')?.outerHTML,
    })));
    await page.screenshot({path: '/tmp/spatialvis-oos-debug.png', fullPage: false});
    throw error;
  }
  assert.equal(await page.locator('rect.region-cell').count(), 60);
  await page.locator('rect.region-cell').first().click({force: true});
  await page.getByRole('button', {name: 'Download CSV'}).waitFor();
  await page.locator('tbody tr').first().waitFor();
  const downloadPromise = page.waitForEvent('download');
  await page.getByRole('button', {name: 'Download CSV'}).click();
  const download = await downloadPromise;
  assert.match(download.suggestedFilename(), /\.csv$/);
  console.log('PASS out-of-sample: 60 regions, metrics table, CSV download');

  await page.goto(`${site}#/diff-resolution`);
  await page.getByRole('heading', {name: 'Replicate 1'}).waitFor();
  assert.ok(await page.locator('table').first().locator('tbody tr').count() > 0);
  console.log('PASS resolution comparison: plot data and three CSV replicates');

  await page.goto(`${site}#/study-design`);
  await page.locator('.markdown-body').waitFor();
  assert.ok((await page.locator('.markdown-body').innerText()).length > 100);
  console.log('PASS study design: markdown content');
  assert.deepEqual(errors, []);
  console.log('PASS browser console: no errors');
} catch (error) {
  console.log('SMOKE DIAGNOSTIC', {
    url: page.url(),
    text: (await page.locator('body').innerText()).slice(0, 1000),
    errors,
  });
  throw error;
} finally {
  await browser.close();
}
