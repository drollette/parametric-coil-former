/**
 * W7HAK Coil Former – Embeddable plugin loader.
 *
 * Drop this script tag into any page on w7hak.com (or anywhere else):
 *
 *   <div id="coil-former"></div>
 *   <script src="https://coil-former.pages.dev/embed.js"
 *           data-target="#coil-former"></script>
 *
 * It creates an iframe pointing at the hosted app, sized to fit.
 */
(function () {
  'use strict';

  const script = document.currentScript;
  const selector = script.getAttribute('data-target') || '#coil-former';
  const width = script.getAttribute('data-width') || '100%';
  const height = script.getAttribute('data-height') || '700px';

  // Determine base URL (same origin as this script)
  const base = new URL('.', script.src).href;

  function mount() {
    const target = document.querySelector(selector);
    if (!target) {
      console.warn('[CoilFormer] Target element not found:', selector);
      return;
    }

    const iframe = document.createElement('iframe');
    iframe.src = base + 'index.html';
    iframe.style.width = width;
    iframe.style.height = height;
    iframe.style.border = 'none';
    iframe.style.borderRadius = '8px';
    iframe.style.overflow = 'hidden';
    iframe.setAttribute('loading', 'lazy');
    iframe.setAttribute('title', 'W7HAK Parametric Coil Former');
    iframe.setAttribute('allow', 'clipboard-write');

    target.appendChild(iframe);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', mount);
  } else {
    mount();
  }
})();
