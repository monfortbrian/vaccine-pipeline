"use client";
import React, { useState, useMemo } from "react";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import type { Epitope } from "@/types";

type SK = "sequence" | "epitope_type" | "hla_allele" | "ic50_nm" | "confidence";

export function EpitopeTable({ epitopes }: { epitopes: Epitope[] }) {
  const [search, setSearch] = useState(""); const [tf, setTf] = useState("All");
  const [sk, setSk] = useState<SK>("ic50_nm"); const [sd, setSd] = useState<"asc" | "desc">("asc");
  const sort = (k: SK) => { if (sk === k) setSd(d => d === "asc" ? "desc" : "asc"); else { setSk(k); setSd("asc"); } };

  const filtered = useMemo(() => {
    let items = [...epitopes];
    if (tf !== "All") items = items.filter(e => e.epitope_type === tf);
    if (search) { const q = search.toLowerCase(); items = items.filter(e => e.sequence.toLowerCase().includes(q) || e.hla_allele.toLowerCase().includes(q)); }
    items.sort((a, b) => { let av: any = a[sk], bv: any = b[sk]; if (av == null) av = sd === "asc" ? Infinity : -Infinity; if (bv == null) bv = sd === "asc" ? Infinity : -Infinity; return typeof av === "string" ? (sd === "asc" ? av.localeCompare(bv) : bv.localeCompare(av)) : (sd === "asc" ? av - bv : bv - av); });
    return items;
  }, [epitopes, search, tf, sk, sd]);

  return (
    <Card><CardHeader><CardTitle>All Epitopes</CardTitle><CardDescription>Sortable and filterable</CardDescription></CardHeader><CardContent>
      <div className="flex flex-wrap items-center gap-3 mb-4">
        <Input placeholder="Search..." value={search} onChange={e => setSearch(e.target.value)} className="max-w-xs h-9 text-sm" />
        <div className="flex gap-1">{["All", "CTL", "HTL", "B-cell"].map(t => <Button key={t} variant={tf === t ? "default" : "outline"} size="sm" onClick={() => setTf(t)}>{t}</Button>)}</div>
        <span className="ml-auto text-xs text-muted-foreground">{filtered.length} epitopes</span>
      </div>
      <div className="overflow-x-auto rounded-md border"><table className="w-full text-sm"><thead><tr className="border-b bg-muted/50">
        {([["sequence", "Sequence"], ["epitope_type", "Type"], ["hla_allele", "HLA"], ["ic50_nm", "IC50 (nM)"], ["confidence", "Confidence"]] as [SK, string][]).map(([k, l]) =>
          <th key={k} onClick={() => sort(k)} className="cursor-pointer select-none px-4 py-3 text-left text-xs font-medium text-muted-foreground hover:text-foreground">{l}{sk === k ? (sd === "asc" ? " ↑" : " ↓") : ""}</th>
        )}
        <th className="px-4 py-3 text-left text-xs font-medium text-muted-foreground">Safety</th>
      </tr></thead><tbody>
          {filtered.length === 0 ? <tr><td colSpan={6} className="py-8 text-center text-muted-foreground">No matches</td></tr> :
            filtered.map((ep, i) => <tr key={`${ep.sequence}-${i}`} className="border-b last:border-0 hover:bg-muted/50 transition-colors">
              <td className="px-4 py-3 font-mono text-xs tracking-wider">{ep.sequence}</td>
              <td className="px-4 py-3"><Badge variant={ep.epitope_type === "CTL" ? "default" : "secondary"}>{ep.epitope_type}</Badge></td>
              <td className="px-4 py-3 font-mono text-xs">{ep.hla_allele}</td>
              <td className="px-4 py-3 font-mono text-xs tabular-nums">{ep.ic50_nm != null ? ep.ic50_nm.toFixed(1) : "-"}</td>
              <td className="px-4 py-3"><Badge variant={ep.confidence === "high" ? "default" : "outline"}>{ep.confidence}</Badge></td>
              <td className="px-4 py-3 text-xs">{ep.allergenicity_safe != null && <Badge variant={ep.allergenicity_safe ? "secondary" : "destructive"} className="mr-1">{ep.allergenicity_safe ? "Safe" : "Allergenic"}</Badge>}{ep.toxicity_safe != null && <Badge variant={ep.toxicity_safe ? "secondary" : "destructive"}>{ep.toxicity_safe ? "Safe" : "Toxic"}</Badge>}</td>
            </tr>)}
        </tbody></table></div>
    </CardContent></Card>
  );
}
