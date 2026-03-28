export type InputType = "pathogen" | "uniprot_id" | "sequence";
export type PipelineStatus = "pending" | "running" | "completed" | "failed";
export type PipelineNode = "N3" | "N4" | "N6" | "N7";
export interface PipelineRunRequest { input_type: InputType; input_value: string; max_proteins?: number; protein_name?: string; }
export interface PipelineRunResponse { run_id: string; status: PipelineStatus; }
export interface PipelineStatusResponse { run_id: string; status: PipelineStatus; current_node: PipelineNode | null; progress: number; message: string; started_at: string | null; completed_at: string | null; }
export interface CoverageDetail { mhc_i_pct: number; mhc_ii_pct: number; combined_pct: number; population_label: string; }
export interface Epitope { sequence: string; epitope_type: "CTL" | "HTL" | "B-cell"; hla_allele: string; ic50_nm: number | null; percentile_rank: number | null; confidence: "high" | "medium" | "low"; allergenicity_safe: boolean | null; toxicity_safe: boolean | null; }
export interface Decision { stage: string; decision: string; reasoning: string; per_population?: Record<string, CoverageDetail>; }
export interface Candidate { protein_id: string; protein_name: string; sequence_length: number; ctl_count: number; ctl_strong: number; htl_count: number; bcell_count: number; global_coverage_pct: number; african_coverage_pct: number; epitopes: Epitope[]; decisions: Decision[]; coverage_detail: Record<string, CoverageDetail>; }
export interface PipelineTiming { total_seconds: number; n3_tcell: number; n4_bcell: number; n6_safety: number; n7_coverage: number; }
export interface PipelineResults { run_id: string; status: PipelineStatus; timing: PipelineTiming; candidates: Candidate[]; }
export interface RunSummary { id: string; pathogen_name: string | null; input_type: InputType; status: PipelineStatus; created_at: string; completed_at: string | null; epitope_count?: number; global_coverage?: number; }
export const PIPELINE_NODES = [
  { id: "N3" as const, label: "T-Cell Epitopes", description: "Predicting CTL & HTL epitopes via IEDB", icon: "🧬" },
  { id: "N4" as const, label: "B-Cell Epitopes", description: "Identifying antibody-binding regions", icon: "🔬" },
  { id: "N6" as const, label: "Safety Screening", description: "Checking allergenicity & toxicity", icon: "🛡️" },
  { id: "N7" as const, label: "Population Coverage", description: "Calculating global HLA coverage", icon: "🌍" },
];
