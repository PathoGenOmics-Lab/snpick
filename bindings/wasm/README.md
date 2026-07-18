# snpick (WebAssembly bindings)

Run snpick's variable-site classification in the browser or Node via WebAssembly.

The `snpick` library is compiled with `default-features = false`, which drops rayon (WASM has no
threads) and memmap (no mmap), using the sequential in-memory scan instead.

## Build

```bash
cargo install wasm-pack
cd bindings/wasm
wasm-pack build --target web
```

## Usage

```js
import init, { variable_site_count, fconst } from "./pkg/snpick_wasm.js";

await init();
const fasta = new TextEncoder().encode(">a\nACGT\n>b\nACGA\n");
console.log(variable_site_count(fasta, false)); // 1
console.log(fconst(fasta, false));              // Uint32Array [1, 1, 1, 0]
```
