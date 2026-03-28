import type {
  PipelineRunRequest, PipelineRunResponse,
  PipelineStatusResponse, PipelineResults, RunSummary,
} from "@/types";

const API = process.env.NEXT_PUBLIC_API_URL || "";

async function f<T>(url: string, opts?: RequestInit): Promise<T> {
  try {
    const r = await fetch(url, { headers: { "Content-Type": "application/json" }, ...opts });
    if (!r.ok) throw new Error(`API returned ${r.status}: ${await r.text()}`);
    return r.json();
  } catch (err) {
    if (err instanceof TypeError && err.message === "Failed to fetch") {
      throw new Error("Backend not reachable. Start FastAPI on port 8001.");
    }
    throw err;
  }
}

export const api = {
  health: () => f<{ status: string }>(`${API}/api/health`),
  startRun: (body: PipelineRunRequest) =>
    f<PipelineRunResponse>(`${API}/api/pipeline/run`, { method: "POST", body: JSON.stringify(body) }),
  getStatus: (id: string) =>
    f<PipelineStatusResponse>(`${API}/api/pipeline/status/${id}`),
  getResults: (id: string) =>
    f<PipelineResults>(`${API}/api/pipeline/results/${id}`),
  listRuns: async () => {
    try {
      const runs = await f<any[]>(`${API}/api/runs`);
      // Normalize: API might return run_id instead of id
      return runs.map((r: any) => ({
        ...r,
        id: r.id || r.run_id,
        pathogen_name: r.pathogen_name || r.input_value || null,
        global_coverage: r.global_coverage ?? r.hla_coverage_global ?? null,
      })) as RunSummary[];
    } catch {
      return [];
    }
  },
};