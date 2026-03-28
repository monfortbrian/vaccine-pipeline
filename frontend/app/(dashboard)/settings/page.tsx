"use client";

import React, { useState } from "react";
import { useRouter } from "next/navigation";
import { useAuth } from "@/components/auth-provider";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Separator } from "@/components/ui/separator";
import { Badge } from "@/components/ui/badge";

export default function SettingsPage() {
  const { user, signOut, updateDisplayName } = useAuth();
  const router = useRouter();
  const [name, setName] = useState(user?.user_metadata?.name || user?.email?.split("@")[0] || "");
  const [copied, setCopied] = useState(false);
  const [saving, setSaving] = useState(false);
  const [saved, setSaved] = useState(false);
  const apiKey = `kz_${(user?.id || "demo").replace(/-/g, "").slice(0, 24)}`;

  const handleSave = async () => {
    setSaving(true);
    try {
      await updateDisplayName(name);
      setSaved(true);
      setTimeout(() => setSaved(false), 2000);
    } catch (e) {
      console.error("Failed to update name:", e);
    } finally {
      setSaving(false);
    }
  };

  return (
    <>
      <div><h2 className="text-2xl font-bold tracking-tight">Settings</h2><p className="text-muted-foreground">Manage your account and API access.</p></div>
      <div className="max-w-2xl space-y-6">
        <Card>
          <CardHeader><CardTitle>Profile</CardTitle><CardDescription>Your account information</CardDescription></CardHeader>
          <CardContent className="space-y-4">
            <div className="grid gap-2">
              <Label htmlFor="display-name">Display name</Label>
              <Input id="display-name" value={name} onChange={(e) => setName(e.target.value)} placeholder="Your name" />
              <p className="text-xs text-muted-foreground">This name appears in the sidebar and audit trails.</p>
            </div>
            <div className="grid gap-2">
              <Label>Email</Label>
              <Input value={user?.email || ""} disabled />
            </div>
            <div className="grid gap-2">
              <Label>Organization</Label>
              <Input value={user?.user_metadata?.organization || "-"} disabled />
            </div>
            <div className="grid gap-2">
              <Label>Role</Label>
              <div className="flex items-center gap-2">
                <Input value={user?.user_metadata?.role || "researcher"} disabled className="flex-1" />
                <Badge>Beta</Badge>
              </div>
            </div>
            <Button size="sm" onClick={handleSave} disabled={saving}>
              {saving ? "Saving..." : saved ? "Saved!" : "Save changes"}
            </Button>
          </CardContent>
        </Card>

        <Card>
          <CardHeader><CardTitle>API Key</CardTitle><CardDescription>For programmatic access (coming soon)</CardDescription></CardHeader>
          <CardContent>
            <div className="flex gap-2">
              <Input value={apiKey} readOnly className="font-mono text-xs" />
              <Button variant="outline" size="sm" className="shrink-0" onClick={() => { navigator.clipboard.writeText(apiKey); setCopied(true); setTimeout(() => setCopied(false), 2000); }}>
                {copied ? "Copied!" : "Copy"}
              </Button>
            </div>
          </CardContent>
        </Card>

        <Separator />
        <div className="flex justify-between items-center">
          <p className="text-sm text-muted-foreground">Need help? <a href="https://kozi-ai.com" target="_blank" rel="noopener noreferrer" className="underline underline-offset-4 hover:text-primary">kozi-ai.com</a></p>
          <Button variant="outline" onClick={async () => { await signOut(); router.push("/login"); }}>Log out</Button>
        </div>
      </div>
    </>
  );
}
