"use client";
import { cn } from "@/lib/utils";
import { PIPELINE_NODES, type PipelineNode } from "@/types";
import { Card, CardContent } from "@/components/ui/card";
import { Progress } from "@/components/ui/progress";
import { Badge } from "@/components/ui/badge";
import { IconDna, IconMicroscope, IconShield, IconWorld } from "@tabler/icons-react";

const ICONS: Record<string, React.ElementType> = { N3: IconDna, N4: IconMicroscope, N6: IconShield, N7: IconWorld };
type S = "waiting" | "running" | "done" | "error";
function st(n: PipelineNode, cur: PipelineNode | null, ps: string): S {
  if (ps === "failed" && n === cur) return "error";
  if (ps === "completed") return "done";
  const o: PipelineNode[] = ["N3","N4","N6","N7"];
  const ci = cur ? o.indexOf(cur) : -1, mi = o.indexOf(n);
  return mi < ci ? "done" : mi === ci ? "running" : "waiting";
}
export function PipelineNodes({ currentNode, status, progress, message }: { currentNode: PipelineNode | null; status: string; progress: number; message: string; }) {
  return (
    <div className="space-y-6">
      <div className="grid grid-cols-1 gap-3 sm:grid-cols-2 lg:grid-cols-4">
        {PIPELINE_NODES.map((node) => { const s = st(node.id, currentNode, status); const Icon = ICONS[node.id] || IconDna; return (
          <Card key={node.id} className={cn("transition-all duration-300", s==="running"&&"ring-2 ring-primary shadow-md", s==="done"&&"opacity-70", s==="waiting"&&"opacity-40", s==="error"&&"ring-2 ring-destructive")}>
            <CardContent className="p-4">
              <div className="flex items-center justify-between mb-3"><div className="rounded-lg bg-muted p-2"><Icon className="size-4 text-muted-foreground" /></div><Badge variant={s==="running"?"default":s==="done"?"secondary":s==="error"?"destructive":"outline"} className="text-[10px]">{s==="running"?"Active":s==="done"?"Done":s==="error"?"Error":"Pending"}</Badge></div>
              <p className="font-medium text-sm">{node.label}</p><p className="text-xs text-muted-foreground mt-1 leading-relaxed">{node.description}</p>
              {s==="running"&&<div className="mt-3 h-1 w-full overflow-hidden rounded-full bg-muted"><div className="h-full bg-primary animate-pulse rounded-full w-2/3"/></div>}
            </CardContent>
          </Card>
        ); })}
      </div>
      <div className="space-y-2">
        <div className="flex items-center justify-between text-sm"><span className="text-muted-foreground">{message||"Initializing..."}</span><span className="font-mono text-xs tabular-nums">{Math.round(progress*100)}%</span></div>
        <Progress value={progress*100} className="h-2" />
      </div>
    </div>
  );
}
