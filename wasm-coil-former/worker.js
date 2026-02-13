/**
 * Cloudflare Worker entry point.
 *
 * Serves the static site from the `dist/` directory.
 * The WASM coil-former runs entirely client-side; the worker simply serves
 * the HTML, JS, and .wasm assets with correct MIME types and caching headers.
 */

const MIME_TYPES = {
  '.html': 'text/html;charset=UTF-8',
  '.js':   'application/javascript',
  '.wasm': 'application/wasm',
  '.css':  'text/css',
  '.json': 'application/json',
  '.png':  'image/png',
  '.svg':  'image/svg+xml',
  '.ico':  'image/x-icon',
};

export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    let path = url.pathname;

    // Default to index.html
    if (path === '/' || path === '') {
      path = '/index.html';
    }

    // CORS headers for plugin embedding on w7hak.com
    const corsHeaders = {
      'Access-Control-Allow-Origin': '*',
      'Access-Control-Allow-Methods': 'GET, OPTIONS',
      'Access-Control-Allow-Headers': 'Content-Type',
    };

    if (request.method === 'OPTIONS') {
      return new Response(null, { status: 204, headers: corsHeaders });
    }

    // Try to serve from the static site bucket
    try {
      const asset = await env.ASSETS.fetch(request);
      if (asset.status === 200) {
        const ext = path.substring(path.lastIndexOf('.'));
        const contentType = MIME_TYPES[ext] || 'application/octet-stream';

        const headers = new Headers(asset.headers);
        headers.set('Content-Type', contentType);
        headers.set('Cache-Control', 'public, max-age=3600');
        Object.entries(corsHeaders).forEach(([k, v]) => headers.set(k, v));

        // .wasm files need special headers for streaming compilation
        if (ext === '.wasm') {
          headers.set('Cache-Control', 'public, max-age=86400');
        }

        return new Response(asset.body, { status: 200, headers });
      }
    } catch (_) {
      // fall through to 404
    }

    return new Response('Not Found', {
      status: 404,
      headers: { 'Content-Type': 'text/plain', ...corsHeaders },
    });
  },
};
