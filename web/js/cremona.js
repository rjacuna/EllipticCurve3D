/* cremona.js -- Cremona's tables of elliptic curves over Q (conductor <= 500,000), served as JSON shards of
 * 10,000 conductors under data/cremona/ (see data/make_cremona.py) and fetched on demand.
 *
 * Browser: the global Cremona.  Node: require() it and point Cremona.fetchJSON at the files.
 *   await Cremona.lookup('11.a1')  ->  { conductor, cremona, lmfdb, ainvs, rank, torsion }  |  null (no such curve)  |  { missing: true, N }
 * rank is the rank of E(Q) and torsion the structure of its torsion subgroup ([2, 4] for Z/2 x Z/4).
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.Cremona = factory();
})(typeof self !== 'undefined' ? self : this, function () {
'use strict';
const C = {
  base: 'data/cremona/',
  index: null,
  shards: new Map(),                                        // file -> shard object
  pending: new Map(),                                       // file -> promise while loading
  fetchJSON: async url => { const r = await fetch(url); if (!r.ok) throw new Error(`HTTP ${r.status} for ${url}`); return r.json(); },
};
C.loadIndex = async function () {
  if (!C.index) C.index = await C.fetchJSON(C.base + 'index.json');
  return C.index;
};
C.maxConductor = () => (C.index ? C.index.max_conductor : null);
// the shard holding conductor N, fetching it if needed; null when N is beyond the tables
C.ensure = async function (N) {
  const idx = await C.loadIndex();
  const sh = idx.shards.find(s => N >= s.lo && N <= s.hi);
  if (!sh) return null;
  if (C.shards.has(sh.file)) return C.shards.get(sh.file);
  if (!C.pending.has(sh.file)) {
    C.pending.set(sh.file, C.fetchJSON(C.base + sh.file).then(d => { C.shards.set(sh.file, d); C.pending.delete(sh.file); return d; },
                                                                 e => { C.pending.delete(sh.file); throw e; }));
  }
  return C.pending.get(sh.file);
};
const LABEL = /^(\d+)(\.?)([a-z]+)(\d+)$/i;
// "11a2" (Cremona) or "11.a1" (LMFDB) -> { N, lmfdb, letter, num }
C.parseLabel = s => { const m = LABEL.exec(String(s).trim()); return m ? { N: Number(m[1]), lmfdb: m[2] === '.', letter: m[3].toLowerCase(), num: Number(m[4]) } : null; };
const rec = (N, cl, cnum, ll, c) => ({ conductor: N, cremona: `${N}${cl}${cnum}`, lmfdb: `${N}.${ll}${c[0]}`, ainvs: c.slice(1, 6),
                                        rank: c[6] === undefined ? null : c[6], torsion: c[7] || [] });
C.lookupIn = function (shard, lab) {
  const classes = shard[String(lab.N)];
  if (!classes) return null;
  for (const [cl, ll, curves] of classes) {
    if (lab.lmfdb ? ll === lab.letter : cl === lab.letter) {
      for (let i = 0; i < curves.length; i++) if (lab.lmfdb ? curves[i][0] === lab.num : i + 1 === lab.num) return rec(lab.N, cl, i + 1, ll, curves[i]);
      return null;
    }
  }
  return null;
};
C.lookup = async function (label) {
  const lab = C.parseLabel(label);
  if (!lab) return null;
  const shard = await C.ensure(lab.N);
  return shard ? C.lookupIn(shard, lab) : { missing: true, N: lab.N };
};
// identify a curve by its a-invariants among the shards already loaded
C.findByAinvs = function (ainvs) {
  const key = ainvs.map(Number).join(',');
  for (const shard of C.shards.values())
    for (const N in shard)
      for (const [cl, ll, curves] of shard[N])
        for (let i = 0; i < curves.length; i++) if (curves[i].slice(1, 6).join(',') === key) return rec(Number(N), cl, i + 1, ll, curves[i]);
  return null;
};
return C;
});
