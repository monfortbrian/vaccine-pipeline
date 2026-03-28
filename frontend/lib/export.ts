import type { Candidate } from "@/types";
function dl(blob: Blob, name: string) { const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(a.href); }
export function downloadJSON(data: unknown, name: string) { dl(new Blob([JSON.stringify(data, null, 2)], { type: "application/json" }), name); }
export function downloadCSV(candidates: Candidate[], name: string) {
  const h = ["protein","protein_id","epitope_sequence","type","hla_allele","ic50_nm","confidence","allergenicity_safe","toxicity_safe"];
  const rows = candidates.flatMap(c => c.epitopes.map(e => [c.protein_name,c.protein_id,e.sequence,e.epitope_type,e.hla_allele,e.ic50_nm!=null?String(e.ic50_nm):"",e.confidence,e.allergenicity_safe!=null?String(e.allergenicity_safe):"",e.toxicity_safe!=null?String(e.toxicity_safe):""]));
  dl(new Blob([h.join(",")+"\n"+rows.map(r=>r.map(c=>`"${c}"`).join(",")).join("\n")], { type: "text/csv" }), name);
}
