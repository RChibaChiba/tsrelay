"""選択 GUI: dock の候補をブラウザで並べ、選んだ ID を返す。

Python 標準の ThreadingHTTPServer を 127.0.0.1 で起動する（WSL2 の localhost は Windows のブラウザから見える）。
"""
import json
import shutil
import subprocess
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

INDEX = r"""<!doctype html>
<html lang="ja"><head><meta charset="utf-8"><title>tsrelay: dock 候補の選択</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.5/3Dmol-min.js"></script>
<style>
body{font-family:system-ui,sans-serif;margin:16px;background:#fafafa;color:#222}
header{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:12px}
button{padding:6px 12px;border:1px solid #888;border-radius:6px;background:#fff;cursor:pointer}
#ok{background:#1d4ed8;color:#fff;border-color:#1d4ed8}
#grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:12px}
.card{background:#fff;border:2px solid #ddd;border-radius:8px;padding:8px}
.card.on{border-color:#1d4ed8}
.view{width:100%;height:260px;position:relative}
.meta{font-size:13px;margin-top:6px;display:flex;justify-content:space-between;align-items:center}
</style></head><body>
<header><b>dock 候補 __N__ 個</b>
<button onclick="sortBy()">エネルギー順</button>
<button onclick="setAll(true)">全選択</button><button onclick="setAll(false)">全解除</button>
<button id="ok" onclick="submit()">確定</button><span id="msg"></span></header>
<div id="grid"></div>
<script>
const C = __DATA__;
const META = __META__;
const grid = document.getElementById("grid");
function card(c){
  const d = document.createElement("div"); d.className = "card"; d.id = "card-" + c.id;
  d.innerHTML = `<div class="view" id="v-${c.id}"></div><div class="meta"><label>
    <input type="checkbox" data-id="${c.id}" onchange="mark(this)"> <b>${c.id}</b></label>
    <span>ΔE ${c.rel_kcal.toFixed(1)} kcal/mol / directed–M ${c.directed_metal.map(x=>x.toFixed(2)).join(", ")} Å</span></div>`;
  grid.appendChild(d);
  const v = $3Dmol.createViewer(document.getElementById("v-" + c.id), {backgroundColor: "white"});
  v.addModel(c.xyz, "xyz");
  v.setStyle({}, {stick: {radius: 0.15}, sphere: {scale: 0.22}});
  v.setStyle({index: META.core}, {stick: {radius: 0.15, color: "#888"}, sphere: {scale: 0.3, color: "#888"}});
  v.setStyle({index: META.directed}, {stick: {radius: 0.15}, sphere: {scale: 0.3, color: "#e11d48"}});
  v.zoomTo(); v.render();
}
function mark(cb){ document.getElementById("card-" + cb.dataset.id).classList.toggle("on", cb.checked); }
function setAll(on){ document.querySelectorAll("input[type=checkbox]").forEach(cb => { cb.checked = on; mark(cb); }); }
function sortBy(){ [...C].sort((a,b) => a.rel_kcal - b.rel_kcal).forEach(c => grid.appendChild(document.getElementById("card-" + c.id))); }
async function submit(){
  const ids = [...document.querySelectorAll("input[type=checkbox]:checked")].map(cb => cb.dataset.id);
  const r = await fetch("/select", {method: "POST", headers: {"Content-Type": "application/json"}, body: JSON.stringify({ids})});
  document.getElementById("msg").textContent = r.ok ? `${ids.length} 個で確定しました。このタブは閉じてかまいません` : "送信に失敗しました";
  document.querySelectorAll("button").forEach(b => b.disabled = true);
}
C.forEach(card);
</script></body></html>
"""


def serve(d, cands, core, directed, port=8000):
    """候補を表示し、確定した ID のリストを返す（空なら何も選ばれなかった）。"""
    d = Path(d)
    data = [dict(c, xyz=(d / f"{c['id']}.xyz").read_text()) for c in cands]
    # 3Dmol の index は 0 始まりの原子番号
    meta = dict(core=list(core), directed=list(directed))
    html = (INDEX.replace("__DATA__", json.dumps(data)).replace("__META__", json.dumps(meta))
            .replace("__N__", str(len(cands)))).encode()
    result = {}

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            if self.path != "/":
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(html)))
            self.end_headers()
            self.wfile.write(html)

        def do_POST(self):
            if self.path != "/select":
                self.send_error(404)
                return
            body = json.loads(self.rfile.read(int(self.headers["Content-Length"])) or b"{}")
            known = {c["id"] for c in cands}
            result["ids"] = [i for i in body.get("ids", []) if i in known]
            self.send_response(200)
            self.end_headers()
            threading.Thread(target=self.server.shutdown, daemon=True).start()

        def log_message(self, *args):
            pass

    try:
        srv = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    except OSError:  # 使用中なら空いているポートを OS に選ばせる
        srv = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    url = f"http://localhost:{srv.server_address[1]}/"
    print(f"候補を選んで確定してください: {url}", flush=True)
    if shutil.which("wslview"):
        subprocess.Popen(["wslview", url], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    srv.serve_forever()
    srv.server_close()
    return result.get("ids", [])
