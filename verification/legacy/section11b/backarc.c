/* Read digraphs in directg -T format ("nv ne a b a b ..."), keep those that are strongly connected, have a
   one-way arc, and pass the back-arc test of Theorem 10.6 (all directed paths with the same vertex set have the
   same number of back arcs). Print survivors as the same text line.
   Used for the pruned cores: geng -c -q 7 | directg -T -q | ./backarc > cores7.txt (107 minutes). */
#include <stdio.h>
#include <string.h>
static int n, adj[16][16], outl[16][16], outc[16];
static int backv[1 << 12];
static int ok;
static int path[16];
static void dfs(int len, int mask, int back) {
    if (!ok) return;
    if (backv[mask] < 0) backv[mask] = back; else if (backv[mask] != back) { ok = 0; return; }
    int v = path[len - 1];
    for (int i = 0; i < outc[v]; i++) {
        int w = outl[v][i];
        if (mask >> w & 1) continue;
        int b = 0;
        for (int j = 0; j < len; j++) if (adj[w][path[j]]) b++;
        path[len] = w;
        dfs(len + 1, mask | (1 << w), back + b);
        if (!ok) return;
    }
}
int main(void) {
    int ne;
    static int ed[400];
    while (scanf("%d %d", &n, &ne) == 2) {
        memset(adj, 0, sizeof adj);
        for (int i = 0; i < ne; i++) { scanf("%d %d", &ed[2*i], &ed[2*i+1]); adj[ed[2*i]][ed[2*i+1]] = 1; }
        /* strong connectivity via Warshall */
        int r[16][16];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) r[i][j] = (i == j) || adj[i][j];
        for (int k = 0; k < n; k++) for (int i = 0; i < n; i++) if (r[i][k]) for (int j = 0; j < n; j++) if (r[k][j]) r[i][j] = 1;
        int sc = 1; for (int i = 0; i < n && sc; i++) for (int j = 0; j < n; j++) if (!r[i][j]) { sc = 0; break; }
        if (!sc) continue;
        int oneway = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (adj[i][j] && !adj[j][i]) oneway = 1;
        if (!oneway) continue;
        for (int i = 0; i < n; i++) { outc[i] = 0; for (int j = 0; j < n; j++) if (adj[i][j]) outl[i][outc[i]++] = j; }
        for (int m = 0; m < (1 << n); m++) backv[m] = -1;
        ok = 1;
        for (int s = 0; s < n && ok; s++) { path[0] = s; dfs(1, 1 << s, 0); }
        if (!ok) continue;
        printf("%d %d", n, ne);
        for (int i = 0; i < ne; i++) printf(" %d %d", ed[2*i], ed[2*i+1]);
        printf("\n");
    }
    return 0;
}
