"use client";

import React, { useState, useEffect, useRef, useCallback } from "react";
import { useRouter } from "next/navigation";
import { useAuth } from "@/components/auth-provider";
import { IconPlayerPlay, IconFlask, IconDna, IconMicroscope, IconShield, IconWorld, IconAlertTriangle } from "@tabler/icons-react";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Textarea } from "@/components/ui/textarea";
import { Switch } from "@/components/ui/switch";
import { Separator } from "@/components/ui/separator";
import { Badge } from "@/components/ui/badge";
import { PipelineNodes } from "@/components/pipeline-nodes";
import { api } from "@/lib/api";
import type { InputType, PipelineStatusResponse } from "@/types";

const INPUT_TYPES: { value: InputType; label: string; desc: string }[] = [
  { value: "pathogen", label: "Pathogen", desc: "Search by name" },
  { value: "uniprot_id", label: "UniProt ID", desc: "Single protein" },
  { value: "sequence", label: "Sequence", desc: "Paste amino acids" },
];

export default function PlaygroundPage() {
  const router = useRouter();
  const { user } = useAuth();
  const firstName = user?.user_metadata?.name?.split(" ")[0] || user?.email?.split("@")[0] || "Researcher";

  const [inputType, setInputType] = useState<InputType>("pathogen");
  const [inputValue, setInputValue] = useState("");
  const [proteinName, setProteinName] = useState("");
  const [maxProteins, setMaxProteins] = useState(3);
  const [safety, setSafety] = useState(true);
  const [coverage, setCoverage] = useState(true);
  const [running, setRunning] = useState(false);
  const [ps, setPs] = useState<PipelineStatusResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const poll = useRef<ReturnType<typeof setInterval> | null>(null);

  const stop = useCallback(() => { if (poll.current) { clearInterval(poll.current); poll.current = null; } }, []);
  const start = useCallback((id: string) => {
    stop();
    poll.current = setInterval(async () => {
      try {
        const s = await api.getStatus(id); setPs(s);
        if (s.status === "completed") { stop(); setRunning(false); router.push(`/results/${id}`); }
        else if (s.status === "failed") { stop(); setRunning(false); setError(s.message || "Pipeline failed"); }
      } catch (e) { console.error(e); }
    }, 2000);
  }, [stop, router]);
  useEffect(() => () => stop(), [stop]);

  const run = async () => {
    if (!inputValue.trim()) return;
    setError(null); setRunning(true); setPs(null);
    try {
      const body = { input_type: inputType, input_value: inputValue.trim() } as any;
      if (inputType === "pathogen") body.max_proteins = maxProteins;
      if (inputType === "sequence" && proteinName) body.protein_name = proteinName;
      const r = await api.startRun(body as Parameters<typeof api.startRun>[0]);
      setPs({ run_id: r.run_id, status: "pending", current_node: null, progress: 0, message: "Pipeline queued...", started_at: null, completed_at: null });
      start(r.run_id);
    } catch (e: unknown) { setError(e instanceof Error ? e.message : "Failed"); setRunning(false); }
  };

  const ph: Record<InputType, string> = { pathogen: "e.g. Mycobacterium tuberculosis", uniprot_id: "e.g. P9WNK7", sequence: "MTEQQWNFAGIEAAASAIQ..." };

  return (
    <>
      <div className="flex items-center justify-between">
        <div>
          <h2 className="text-2xl font-bold tracking-tight">Good day, {firstName} 👋</h2>
          <p className="text-muted-foreground mt-1">Analyze any pathogen, protein, or amino acid sequence.</p>
        </div>
        <Badge variant="outline">Beta</Badge>
      </div>

      <div className="grid gap-6 lg:grid-cols-[400px_1fr]">
        <div className="flex flex-col gap-4">
          <Card>
            <CardHeader><CardTitle>Input Configuration</CardTitle><CardDescription>Select input type and set parameters</CardDescription></CardHeader>
            <CardContent className="space-y-6">
              <div className="grid grid-cols-3 gap-2">
                {INPUT_TYPES.map((t) => (
                  <button key={t.value} type="button" onClick={() => { setInputType(t.value); setInputValue(""); }}
                    className={`flex flex-col items-center gap-1 rounded-lg border-2 px-3 py-3 text-center transition-all cursor-pointer ${inputType === t.value ? "border-primary bg-primary/5 text-primary" : "border-muted bg-transparent text-muted-foreground hover:bg-accent hover:text-accent-foreground"}`}>
                    <span className="text-xs font-semibold">{t.label}</span>
                    <span className="text-[10px] leading-tight opacity-70">{t.desc}</span>
                  </button>
                ))}
              </div>
              <Separator />
              {inputType === "sequence" ? (
                <div className="space-y-4">
                  <div className="grid gap-2"><Label>Amino acid sequence</Label><Textarea className="font-mono text-xs min-h-[120px] resize-none" placeholder={ph.sequence} value={inputValue} onChange={(e) => setInputValue(e.target.value)} disabled={running} /></div>
                  <div className="grid gap-2"><Label>Protein name <span className="text-muted-foreground font-normal">(optional)</span></Label><Input placeholder="Custom protein" value={proteinName} onChange={(e) => setProteinName(e.target.value)} disabled={running} /></div>
                </div>
              ) : (
                <div className="grid gap-2"><Label>{inputType === "pathogen" ? "Pathogen name" : "UniProt ID"}</Label><Input placeholder={ph[inputType]} value={inputValue} onChange={(e) => setInputValue(e.target.value)} disabled={running} className={inputType === "uniprot_id" ? "font-mono" : ""} /></div>
              )}
              {inputType === "pathogen" && (<><Separator /><div className="space-y-3"><div className="flex items-center justify-between"><Label>Max proteins to analyze</Label><span className="rounded-md bg-primary text-primary-foreground px-2.5 py-0.5 font-mono text-xs font-bold">{maxProteins}</span></div><input type="range" min={1} max={10} step={1} value={maxProteins} onChange={(e) => setMaxProteins(Number(e.target.value))} disabled={running} className="w-full accent-primary h-2 rounded-full cursor-pointer" /><div className="flex justify-between text-[10px] text-muted-foreground"><span>1</span><span>5</span><span>10</span></div></div></>)}
              <Separator />
              <div className="space-y-4">
                <div className="flex items-center justify-between"><div><Label className="text-sm">Safety screening</Label><p className="text-[11px] text-muted-foreground">Allergenicity & toxicity checks</p></div><Switch checked={safety} onCheckedChange={setSafety} disabled={running} /></div>
                <div className="flex items-center justify-between"><div><Label className="text-sm">Population coverage</Label><p className="text-[11px] text-muted-foreground">HLA coverage across 7 populations</p></div><Switch checked={coverage} onCheckedChange={setCoverage} disabled={running} /></div>
              </div>
            </CardContent>
          </Card>
          <Button className="w-full" size="lg" onClick={run} disabled={running || !inputValue.trim()}>
            {running ? <><span className="h-4 w-4 animate-spin rounded-full border-2 border-current border-t-transparent mr-2" />Running pipeline...</> : <><IconPlayerPlay className="size-4 mr-2" />Run Analysis</>}
          </Button>
        </div>
        <div className="flex flex-col gap-4">
          {error && (
            <Card className="border-destructive/50 bg-destructive/5"><CardContent className="flex items-start gap-3 pt-6">
              <div className="rounded-full bg-destructive/10 p-2 shrink-0"><IconAlertTriangle className="size-4 text-destructive" /></div>
              <div><p className="text-sm font-semibold text-destructive">Connection Error</p><p className="text-sm text-muted-foreground mt-1">{error}</p><p className="text-xs text-muted-foreground mt-3">Start your backend:</p><code className="mt-1 block rounded bg-muted px-3 py-2 text-xs font-mono">cd api && uvicorn main:app --port 8001 --reload</code></div>
            </CardContent></Card>
          )}
          {!running && !ps && !error && (
            <div className="rounded-xl border border-dashed flex flex-col items-center justify-center py-24 px-8 text-center">
              <div className="rounded-full bg-muted p-4 mb-6"><IconFlask className="size-8 text-muted-foreground" /></div>
              <h3 className="text-lg font-semibold">Ready to analyze</h3>
              <p className="text-sm text-muted-foreground mt-2 max-w-sm">Enter a pathogen, UniProt ID, or sequence. The pipeline predicts epitopes, screens for safety, and calculates population coverage.</p>
              <div className="grid grid-cols-4 gap-8 mt-10">
                {[{ icon: IconDna, label: "T-Cell" }, { icon: IconMicroscope, label: "B-Cell" }, { icon: IconShield, label: "Safety" }, { icon: IconWorld, label: "Coverage" }].map((s) => (
                  <div key={s.label} className="flex flex-col items-center gap-2"><div className="rounded-lg bg-muted p-2.5"><s.icon className="size-5 text-muted-foreground" /></div><span className="text-xs text-muted-foreground font-medium">{s.label}</span></div>
                ))}
              </div>
            </div>
          )}
          {ps && (<Card><CardHeader><div className="flex items-center justify-between"><CardTitle>Pipeline Progress</CardTitle><Badge variant="outline" className="font-mono text-xs">{ps.run_id.slice(0, 8)}</Badge></div></CardHeader><CardContent><PipelineNodes currentNode={ps.current_node} status={ps.status} progress={ps.progress} message={ps.message} /></CardContent></Card>)}
        </div>
      </div>
    </>
  );
}