"use client";

import React, { useEffect, useState } from "react";
import { useRouter } from "next/navigation";
import { IconPlus } from "@tabler/icons-react";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Card, CardContent } from "@/components/ui/card";
import { api } from "@/lib/api";
import type { RunSummary } from "@/types";

function formatDate(iso: string): string {
  try {
    const d = new Date(iso);
    if (isNaN(d.getTime())) return "-";
    return d.toLocaleDateString("en-US", {
      month: "short", day: "numeric", year: "numeric",
      hour: "2-digit", minute: "2-digit",
    });
  } catch {
    return "-";
  }
}

export default function HistoryPage() {
  const router = useRouter();
  const [runs, setRuns] = useState<RunSummary[]>([]);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    api.listRuns().then(setRuns).catch(() => setRuns([])).finally(() => setLoading(false));
  }, []);

  return (
    <>
      <div className="flex items-center justify-between">
        <div>
          <h2 className="text-2xl font-bold tracking-tight">Run History</h2>
          <p className="text-muted-foreground">All pipeline runs and their results.</p>
        </div>
        <Button onClick={() => router.push("/playground")}>
          <IconPlus className="mr-2 size-4" />New Run
        </Button>
      </div>

      {loading ? (
        <div className="flex justify-center py-12">
          <div className="h-6 w-6 animate-spin rounded-full border-2 border-primary border-t-transparent" />
        </div>
      ) : runs.length === 0 ? (
        <Card>
          <CardContent className="flex flex-col items-center py-20">
            <div className="rounded-full bg-muted p-4 mb-4">
              <IconPlus className="size-6 text-muted-foreground" />
            </div>
            <p className="font-medium">No runs yet</p>
            <p className="text-sm text-muted-foreground mt-1">Start your first analysis in the Playground.</p>
            <Button className="mt-6" onClick={() => router.push("/playground")}>Go to Playground</Button>
          </CardContent>
        </Card>
      ) : (
        <div className="rounded-lg border overflow-hidden">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b bg-muted/50">
                <th className="text-left px-4 py-3 text-xs font-medium text-muted-foreground">Date</th>
                <th className="text-left px-4 py-3 text-xs font-medium text-muted-foreground">Pathogen</th>
                <th className="text-left px-4 py-3 text-xs font-medium text-muted-foreground">Status</th>
                <th className="text-right px-4 py-3 text-xs font-medium text-muted-foreground">Coverage</th>
                <th className="px-4 py-3"></th>
              </tr>
            </thead>
            <tbody>
              {runs.map((r, index) => (
                <tr
                  key={r.id || `run-${index}`}
                  className="border-b last:border-0 hover:bg-muted/50 cursor-pointer transition-colors"
                  onClick={() => r.status === "completed" && router.push(`/results/${r.id}`)}
                >
                  <td className="px-4 py-3 text-xs text-muted-foreground">{formatDate(r.created_at)}</td>
                  <td className="px-4 py-3 font-medium">{r.pathogen_name || "-"}</td>
                  <td className="px-4 py-3">
                    <Badge variant={r.status === "completed" ? "default" : r.status === "failed" ? "destructive" : "secondary"}>
                      {r.status}
                    </Badge>
                  </td>
                  <td className="px-4 py-3 text-right font-mono text-xs">
                    {r.global_coverage != null ? `${r.global_coverage.toFixed(1)}%` : "-"}
                  </td>
                  <td className="px-4 py-3 text-right">
                    {r.status === "completed" && <Button variant="ghost" size="sm">View →</Button>}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </>
  );
}
