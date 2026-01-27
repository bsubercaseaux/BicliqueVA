import React, { useState, useMemo } from "react";

// ---------- Types ----------

type Clause = number[]; // list of integers, DIMACS-style literal ids

interface CnfFormula {
    numVars: number;
    clauses: Clause[];
}

interface Graph {
    vertices: number[];
    adjacency: Record<number, Set<number>>;
}

interface BvaStepInfo {
    stepIndex: number;
    L: number[];
    R: number[];
    quality: number;
    newVar: number;
    removedEdges: number;
    addedEdges: number;
}

// NEW: priority queue item (vertex + degree)
interface PriorityItem {
    v: number;
    degree: number;
}

// ---------- Utility helpers ----------

function cloneGraph(g: Graph): Graph {
    const adjacency: Record<number, Set<number>> = {};
    for (const [k, v] of Object.entries(g.adjacency)) {
        adjacency[+k] = new Set(v);
    }
    return {
        vertices: [...g.vertices],
        adjacency,
    };
}

function graphDegree(g: Graph, v: number): number {
    return g.adjacency[v]?.size ?? 0;
}

function graphAddVertex(g: Graph, v: number) {
    if (!(v in g.adjacency)) {
        g.adjacency[v] = new Set();
        if (!g.vertices.includes(v)) {
            g.vertices.push(v);
        }
    }
}

function graphAddEdge(g: Graph, u: number, v: number) {
    if (u === v) return;
    graphAddVertex(g, u);
    graphAddVertex(g, v);
    if (!g.adjacency[u].has(v)) {
        g.adjacency[u].add(v);
        g.adjacency[v].add(u);
    }
}

function graphRemoveEdge(g: Graph, u: number, v: number) {
    if (g.adjacency[u]?.has(v)) {
        g.adjacency[u].delete(v);
        g.adjacency[v].delete(u);
    }
}

function countEdges(g: Graph): number {
    let total = 0;
    for (const v of g.vertices) {
        total += g.adjacency[v]?.size ?? 0;
    }
    return total / 2;
}

function bvaQuality(k: number, s: number): number {
    // Q = |L||R| - |L| - |R|
    return k * s - k - s;
}

// ---------- DIMACS parser & graph builder ----------

function parseDimacsCnf(text: string): CnfFormula {
    const lines = text
        .split(/\r?\n/)
        .map((l) => l.trim())
        .filter((l) => l.length > 0 && !l.startsWith("c"));

    let numVars = 0;
    const clauses: Clause[] = [];

    for (const line of lines) {
        if (line.startsWith("p")) {
            const parts = line.split(/\s+/);
            if (parts.length >= 4 && parts[1] === "cnf") {
                numVars = parseInt(parts[2], 10);
            }
        } else {
            const lits = line
                .split(/\s+/)
                .map((x) => parseInt(x, 10))
                .filter((x) => x !== 0);
            if (lits.length > 0) {
                clauses.push(lits);
            }
        }
    }
    if (numVars === 0) {
        for (const cl of clauses) {
            for (const lit of cl) {
                numVars = Math.max(numVars, Math.abs(lit));
            }
        }
    }
    return { numVars, clauses };
}

/**
 * Build a graph from monotone negative 2-clauses:
 * each clause (-u -v) becomes an undirected edge {u,v}.
 */
function buildGraphFromCnf(formula: CnfFormula): Graph {
    const adjacency: Record<number, Set<number>> = {};
    const vertices = new Set<number>();

    for (const clause of formula.clauses) {
        if (clause.length === 2 && clause[0] < 0 && clause[1] < 0) {
            const u = -clause[0];
            const v = -clause[1];
            vertices.add(u);
            vertices.add(v);
            if (!adjacency[u]) adjacency[u] = new Set();
            if (!adjacency[v]) adjacency[v] = new Set();
            adjacency[u].add(v);
            adjacency[v].add(u);
        }
    }

    for (let i = 1; i <= formula.numVars; i++) {
        vertices.add(i);
        if (!adjacency[i]) adjacency[i] = new Set();
    }

    return {
        vertices: Array.from(vertices).sort((a, b) => a - b),
        adjacency,
    };
}

// ---------- Priority queue helpers (by degree) ----------

function buildPriorityQueue(g: Graph): PriorityItem[] {
    const items: PriorityItem[] = [];
    for (const v of g.vertices) {
        const d = graphDegree(g, v);
        if (d > 0) {
            items.push({ v, degree: d });
        }
    }
    // sort: highest degree first, then smaller variable id
    items.sort((a, b) => {
        if (b.degree !== a.degree) return b.degree - a.degree;
        return a.v - b.v;
    });
    return items;
}

// ---------- BVA step (graph-based S-BVA heuristic with PQ) ----------

/**
 * Compute one BVA step as a biclique (L,R) in the current graph.
 * If seedHint is provided, we try to use that as seed (top of the priority queue),
 * falling back to the maximum-degree vertex if that seed is now isolated.
 */
function computeBvaStep(
    graph: Graph,
    seedHint?: number
): { L: number[]; R: number[]; quality: number } | null {
    const degrees = new Map<number, number>();
    const nonIsolated: number[] = [];

    for (const v of graph.vertices) {
        const d = graphDegree(graph, v);
        degrees.set(v, d);
        if (d > 0) nonIsolated.push(v);
    }

    if (nonIsolated.length === 0) {
        return null;
    }

    let seed: number;
    // Try to use the head of the priority queue as in SimpleBVA (Fig. 1, line 3) :contentReference[oaicite:1]{index=1}
    if (seedHint !== undefined && (degrees.get(seedHint) ?? 0) > 0) {
        seed = seedHint;
    } else {
        // fallback: recompute max-degree vertex if head is stale
        seed = nonIsolated[0];
        for (const v of nonIsolated) {
            if ((degrees.get(v) ?? 0) > (degrees.get(seed) ?? 0)) {
                seed = v;
            }
        }
    }

    let L = new Set<number>([seed]);
    let R = new Set<number>(graph.adjacency[seed] ?? new Set());
    let Q = bvaQuality(L.size, R.size);

    while (true) {
        let bestV: number | null = null;
        let bestR: Set<number> | null = null;
        let bestQ = Q;

        for (const v of graph.vertices) {
            if (L.has(v)) continue;
            const neigh = graph.adjacency[v] ?? new Set<number>();
            const newR = new Set<number>();
            for (const u of R) {
                if (neigh.has(u)) newR.add(u);
            }
            if (newR.size === 0) continue;

            const newQ = bvaQuality(L.size + 1, newR.size);
            if (newQ > bestQ) {
                bestQ = newQ;
                bestR = newR;
                bestV = v;
            }
        }

        if (bestV === null) break;
        L.add(bestV);
        R = bestR!;
        Q = bestQ;
    }

    if (L.size >= 2 && Q > 0 && R.size > 0) {
        return {
            L: Array.from(L),
            R: Array.from(R),
            quality: Q,
        };
    }
    return null;
}

/**
 * Apply a BVA step to both graph and CNF.
 */
function applyBvaStep(
    graph: Graph,
    formula: CnfFormula,
    step: { L: number[]; R: number[]; quality: number }
): { newGraph: Graph; newFormula: CnfFormula; info: BvaStepInfo } {
    const { L, R, quality } = step;
    const g = cloneGraph(graph);
    const numVarsBefore = formula.numVars;
    const newVar = numVarsBefore + 1;

    graphAddVertex(g, newVar);

    let removedEdges = 0;
    let addedEdges = 0;

    for (const u of L) {
        if (!g.adjacency[newVar].has(u)) {
            graphAddEdge(g, u, newVar);
            addedEdges++;
        }
    }
    for (const v of R) {
        if (!g.adjacency[newVar].has(v)) {
            graphAddEdge(g, v, newVar);
            addedEdges++;
        }
    }

    for (const u of L) {
        for (const v of R) {
            if (g.adjacency[u]?.has(v)) {
                graphRemoveEdge(g, u, v);
                removedEdges++;
            }
        }
    }

    const newClauses: Clause[] = [];
    const Lset = new Set(L);
    const Rset = new Set(R);

    for (const clause of formula.clauses) {
        if (clause.length === 2 && clause[0] < 0 && clause[1] < 0) {
            const a = -clause[0];
            const b = -clause[1];
            const isInBiclique =
                (Lset.has(a) && Rset.has(b)) || (Lset.has(b) && Rset.has(a));
            if (isInBiclique) {
                continue;
            }
        }
        newClauses.push(clause);
    }

    const yLit = -newVar;
    for (const u of L) {
        newClauses.push([-u, yLit]);
    }
    for (const v of R) {
        newClauses.push([yLit, -v]);
    }

    const newFormula: CnfFormula = {
        numVars: newVar,
        clauses: newClauses,
    };

    const info: BvaStepInfo = {
        stepIndex: -1,
        L: [...L].sort((a, b) => a - b),
        R: [...R].sort((a, b) => a - b),
        quality,
        newVar,
        removedEdges,
        addedEdges,
    };

    return { newGraph: g, newFormula, info };
}

// ---------- Graph visualization ----------

interface GraphViewProps {
    graph: Graph | null;
    L: Set<number>;
    R: Set<number>;
    baseVarCount: number;
}

const GraphView: React.FC<GraphViewProps> = ({ graph, L, R, baseVarCount }) => {
    const width = 640;
    const height = 480;
    const radius = 200;
    const centerX = width / 2;
    const centerY = height / 2;

    if (!graph) {
        return <div className="panel">No graph to display yet.</div>;
    }

    const verts = graph.vertices.slice().sort((a, b) => a - b);
    const positions = new Map<number, { x: number; y: number }>();
    const n = verts.length;

    verts.forEach((v, idx) => {
        const angle = (2 * Math.PI * idx) / n;
        const x = centerX + radius * Math.cos(angle);
        const y = centerY + radius * Math.sin(angle);
        positions.set(v, { x, y });
    });

    const edges: [number, number][] = [];
    for (const u of verts) {
        for (const v of graph.adjacency[u]) {
            if (u < v) edges.push([u, v]);
        }
    }

    const nodeColor = (v: number) => {
        if (v > baseVarCount) {
            return "#2e7d32"; // aux variables -> green
        }
        if (L.has(v)) return "#1976d2"; // L -> blue
        if (R.has(v)) return "#ef6c00"; // R -> orange
        return "#9e9e9e"; // other base vars -> grey
    };

    return (
        <div className="panel">
            <h3>Graph view (monotone 2-CNF part)</h3>
            <svg width={width} height={height} style={{ border: "1px solid #ccc" }}>
                {edges.map(([u, v], idx) => {
                    const pu = positions.get(u)!;
                    const pv = positions.get(v)!;
                    return (
                        <line
                            key={idx}
                            x1={pu.x}
                            y1={pu.y}
                            x2={pv.x}
                            y2={pv.y}
                            stroke="#ddd"
                            strokeWidth={1}
                        />
                    );
                })}
                {verts.map((v) => {
                    const p = positions.get(v)!;
                    return (
                        <g key={v}>
                            <circle
                                cx={p.x}
                                cy={p.y}
                                r={14}
                                fill={nodeColor(v)}
                                stroke="#333"
                                strokeWidth={1}
                            />
                            <text
                                x={p.x}
                                y={p.y + 4}
                                textAnchor="middle"
                                fontSize="10"
                                fill="#fff"
                            >
                                {v}
                            </text>
                        </g>
                    );
                })}
            </svg>
            <p style={{ fontSize: "0.85rem", marginTop: 8 }}>
                Colors: blue = L, orange = R, green = BVA aux vars, grey = other vars.
            </p>
        </div>
    );
};

// ---------- Main App ----------

const App: React.FC = () => {
    const [rawText, setRawText] = useState<string>("");
    const [formula, setFormula] = useState<CnfFormula | null>(null);
    const [baseFormula, setBaseFormula] = useState<CnfFormula | null>(null);
    const [graph, setGraph] = useState<Graph | null>(null);
    const [baseGraph, setBaseGraph] = useState<Graph | null>(null);
    const [steps, setSteps] = useState<BvaStepInfo[]>([]);
    const [status, setStatus] = useState<string>("Load a CNF file to begin.");

    // NEW: priority queue state
    const [priorityQueue, setPriorityQueue] = useState<PriorityItem[]>([]);

    const currentStep = steps.length > 0 ? steps[steps.length - 1] : null;
    const currentL = useMemo(
        () => new Set(currentStep ? currentStep.L : []),
        [currentStep]
    );
    const currentR = useMemo(
        () => new Set(currentStep ? currentStep.R : []),
        [currentStep]
    );

    const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0];
        if (!file) return;
        const reader = new FileReader();
        reader.onload = () => {
            const text = String(reader.result ?? "");
            setRawText(text);
            const parsed = parseDimacsCnf(text);
            const g = buildGraphFromCnf(parsed);
            setFormula(parsed);
            setBaseFormula(parsed);
            setGraph(g);
            setBaseGraph(g);
            setSteps([]);
            setPriorityQueue(buildPriorityQueue(g)); // NEW
            setStatus(
                `Loaded CNF with ${parsed.numVars} vars and ${parsed.clauses.length} clauses.`
            );
        };
        reader.readAsText(file);
    };

    const handleReset = () => {
        if (!baseFormula || !baseGraph) return;
        const g = cloneGraph(baseGraph);
        setFormula({ ...baseFormula, clauses: [...baseFormula.clauses] });
        setGraph(g);
        setSteps([]);
        setPriorityQueue(buildPriorityQueue(g)); // NEW
        setStatus("Reset to original CNF and graph.");
    };

    const runNextBvaStep = () => {
        if (!graph || !formula) {
            setStatus("Please load a CNF file first.");
            return;
        }
        const seedHint = priorityQueue.length > 0 ? priorityQueue[0].v : undefined;
        const step = computeBvaStep(graph, seedHint);
        if (!step) {
            setStatus(
                "No more profitable BVA steps found (no biclique (L,R) with Q>0)."
            );
            return;
        }

        const { newGraph, newFormula, info } = applyBvaStep(graph, formula, step);
        const newStepIndex = steps.length;
        info.stepIndex = newStepIndex;

        setGraph(newGraph);
        setFormula(newFormula);
        setSteps([...steps, info]);
        setPriorityQueue(buildPriorityQueue(newGraph)); // NEW

        setStatus(
            `Applied BVA step #${newStepIndex}: seed=${seedHint ?? "(max-degree)"} introduced y=${info.newVar}, removed ${info.removedEdges} edges, added ${info.addedEdges} edges.`
        );
    };

    const currentEdgeCount = graph ? countEdges(graph) : 0;
    const baseEdgeCount = baseGraph ? countEdges(baseGraph) : 0;
    const baseVarCount = baseFormula?.numVars ?? 0;
    const currentVarCount = formula?.numVars ?? 0;
    const currentClauseCount = formula?.clauses.length ?? 0;
    const baseClauseCount = baseFormula?.clauses.length ?? 0;

    return (
        <div className="app-root">
            <header className="app-header">
                <h1>BVA Visualizer</h1>
                <p>
                    Load a (DIMACS) CNF file, then step through the graph-based BVA
                    heuristic on the monotone 2-CNF part.
                </p>
            </header>

            <main className="app-main">
                <section className="left-column">
                    <div className="panel">
                        <h3>1. Load CNF</h3>
                        <input type="file" accept=".cnf,.dimacs,.txt" onChange={handleFileChange} />
                        <button
                            onClick={handleReset}
                            disabled={!baseFormula || !baseGraph}
                            style={{ marginLeft: 8 }}
                        >
                            Reset
                        </button>
                        <p className="status-text">{status}</p>
                        {formula && (
                            <div className="stats">
                                <h4>Formula stats</h4>
                                <ul>
                                    <li>
                                        Vars: {currentVarCount}{" "}
                                        {baseVarCount > 0 && (
                                            <span className="muted">
                                                (initial {baseVarCount}, added{" "}
                                                {currentVarCount - baseVarCount})
                                            </span>
                                        )}
                                    </li>
                                    <li>
                                        Clauses: {currentClauseCount}{" "}
                                        {baseClauseCount > 0 && (
                                            <span className="muted">
                                                (initial {baseClauseCount}, Δ{" "}
                                                {currentClauseCount - baseClauseCount})
                                            </span>
                                        )}
                                    </li>
                                    <li>
                                        Edges in monotone 2-CNF graph: {currentEdgeCount}{" "}
                                        {baseGraph && (
                                            <span className="muted">
                                                (initial {baseEdgeCount}, Δ{" "}
                                                {currentEdgeCount - baseEdgeCount})
                                            </span>
                                        )}
                                    </li>
                                </ul>
                            </div>
                        )}
                    </div>

                    <div className="panel">
                        <h3>2. BVA control</h3>
                        <button onClick={runNextBvaStep} disabled={!graph || !formula}>
                            Next BVA step
                        </button>
                        {currentStep && (
                            <div style={{ marginTop: 16 }}>
                                <h4>Current step #{currentStep.stepIndex}</h4>
                                <p>
                                    Introduced new variable <strong>y = {currentStep.newVar}</strong>
                                    .
                                </p>
                                <p>
                                    |L| = <strong>{currentStep.L.length}</strong>, |R| ={" "}
                                    <strong>{currentStep.R.length}</strong>, reduction Q(L,R) ={" "}
                                    <strong>{currentStep.quality}</strong>
                                </p>
                                <p>
                                    Edges removed: <strong>{currentStep.removedEdges}</strong>, edges
                                    added: <strong>{currentStep.addedEdges}</strong>.
                                </p>
                                <details style={{ marginTop: 8 }}>
                                    <summary>Show L and R</summary>
                                    <p>
                                        L = {"{"}
                                        {currentStep.L.join(", ")}
                                        {"}"}
                                    </p>
                                    <p>
                                        R = {"{"}
                                        {currentStep.R.join(", ")}
                                        {"}"}
                                    </p>
                                </details>
                            </div>
                        )}
                        {!currentStep && (
                            <p className="muted" style={{ marginTop: 8 }}>
                                No BVA step applied yet.
                            </p>
                        )}
                    </div>

                    {/* NEW: Priority queue view */}
                    <div className="panel">
                        <h3>3. Priority queue (by degree)</h3>
                        {priorityQueue.length === 0 ? (
                            <p className="muted">Queue is empty (no non-isolated vertices).</p>
                        ) : (
                            <div style={{ fontSize: "0.85rem" }}>
                                <p>
                                    Top of queue is{" "}
                                    <strong>{priorityQueue[0].v}</strong> (deg{" "}
                                    {priorityQueue[0].degree}).
                                </p>
                                <div
                                    style={{
                                        maxHeight: 140,
                                        overflow: "auto",
                                        border: "1px solid #eee",
                                        padding: "4px 6px",
                                        borderRadius: 4,
                                        background: "#fafafa",
                                    }}
                                >
                                    {priorityQueue.map((item, idx) => (
                                        <div
                                            key={item.v}
                                            style={{
                                                display: "flex",
                                                justifyContent: "space-between",
                                                padding: "1px 0",
                                                fontWeight: idx === 0 ? 600 : 400,
                                                color: idx === 0 ? "#1565c0" : "#333",
                                            }}
                                        >
                                            <span>
                                                {idx === 0 && "▶ "}v={item.v}
                                            </span>
                                            <span>deg={item.degree}</span>
                                        </div>
                                    ))}
                                </div>
                            </div>
                        )}
                    </div>
                </section>

                <section className="right-column">
                    <GraphView
                        graph={graph}
                        L={currentL}
                        R={currentR}
                        baseVarCount={baseVarCount}
                    />
                </section>
            </main>

            <footer className="app-footer">
                <p>
                    Note: This visualizer uses the monotone negative 2-CNF part of your
                    formula and a graph-based SimpleBVA step (priority queue of literals
                    by |F<sub>ℓ</sub>|, i.e. degree). Other clauses remain in the CNF but
                    are not shown in the graph.
                </p>
            </footer>
        </div>
    );
};

export default App;
