"use client";
import { Bar, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, ReferenceLine, Line, ComposedChart } from "recharts";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import type { CoverageDetail } from "@/types";

const ORDER = ["global", "african", "east_african", "european", "east_asian", "south_asian", "americas"];

export function CoverageChart({ coverageDetail }: { coverageDetail: Record<string, CoverageDetail> }) {
  const data = ORDER.filter(k => coverageDetail[k]).map(k => { const d = coverageDetail[k]; return { name: d.population_label, "MHC-I": +d.mhc_i_pct.toFixed(1), "MHC-II": +d.mhc_ii_pct.toFixed(1), Combined: +d.combined_pct.toFixed(1) }; });
  return (
    <Card>
      <CardHeader><CardTitle>Population Coverage</CardTitle><CardDescription>HLA coverage by population - 80% target</CardDescription></CardHeader>
      <CardContent><div className="h-[360px] w-full"><ResponsiveContainer width="100%" height="100%">
        <ComposedChart data={data} margin={{ top: 8, right: 12, left: 0, bottom: 0 }}>
          <CartesianGrid strokeDasharray="3 3" className="stroke-muted" />
          <XAxis dataKey="name" tick={{ fontSize: 11 }} angle={-20} textAnchor="end" height={60} />
          <YAxis domain={[0, 100]} tick={{ fontSize: 11 }} tickFormatter={v => `${v}%`} />
          <Tooltip formatter={(v: number) => [`${v}%`]} />
          <Legend wrapperStyle={{ fontSize: 12 }} />
          <ReferenceLine y={80} stroke="var(--color-muted-foreground)" strokeDasharray="6 4" label={{ value: "80% target", position: "right", fontSize: 10 }} />
          <Bar dataKey="MHC-I" fill="var(--color-primary)" radius={[2, 2, 0, 0]} barSize={20} />
          <Bar dataKey="MHC-II" fill="var(--color-muted-foreground)" radius={[2, 2, 0, 0]} barSize={20} />
          <Line type="monotone" dataKey="Combined" stroke="var(--color-chart-1)" strokeWidth={2} dot={{ r: 4 }} />
        </ComposedChart>
      </ResponsiveContainer></div></CardContent>
    </Card>
  );
}
