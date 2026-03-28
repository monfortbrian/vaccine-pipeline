"use client";

import React, { useState } from "react";
import { useRouter } from "next/navigation";
import { IconFlask, IconEye, IconEyeOff } from "@tabler/icons-react";
import { cn } from "@/lib/utils";
import { Button } from "@/components/ui/button";
import { Card, CardContent } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { useAuth } from "@/components/auth-provider";
import { AccessDeniedModal } from "@/components/access-denied-modal";

export default function LoginPage() {
  const router = useRouter();
  const { signIn, user } = useAuth();
  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [showPassword, setShowPassword] = useState(false);
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(false);
  const [showAccessDenied, setShowAccessDenied] = useState(false);

  React.useEffect(() => {
    if (user) router.push("/playground");
  }, [user, router]);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!email || !password) { setError("Please fill in all fields."); return; }
    setError("");
    setLoading(true);
    try {
      await signIn(email, password);
      router.push("/playground");
    } catch (err: unknown) {
      const msg = err instanceof Error ? err.message : "Sign in failed";
      if (msg === "ACCESS_DENIED") {
        setShowAccessDenied(true);
      } else if (msg.includes("Invalid login credentials")) {
        setError("Incorrect email or password. Please try again.");
      } else {
        setError(msg);
      }
    } finally {
      setLoading(false);
    }
  };

  return (
    <>
      {showAccessDenied && <AccessDeniedModal onClose={() => setShowAccessDenied(false)} />}
      <div className="flex min-h-svh flex-col items-center justify-center bg-mist-800 p-10 md:p-10">
        <div className="w-full max-w-sm md:max-w-4xl">
          <div className={cn("flex flex-col gap-6")}>
            <Card className="overflow-hidden p-0">
              <CardContent className="grid p-0 md:grid-cols-2">
                <form className="p-6 md:p-8" onSubmit={handleSubmit}>
                  <div className="flex flex-col gap-6">
                    <div className="flex flex-col items-center gap-2 text-center">
                      <div className="flex size-11 items-center justify-center rounded-xl bg-primary text-primary-foreground">
                        <IconFlask className="size-6" />
                      </div>
                      <h1 className="text-2xl font-bold">Welcome back</h1>
                      <p className="text-balance text-sm text-muted-foreground">Sign in to the Kozi AI platform</p>
                    </div>
                    {error && <div className="rounded-lg border border-destructive/50 bg-destructive/10 px-4 py-3 text-sm text-destructive">{error}</div>}
                    <div className="grid gap-4">
                      <div className="grid gap-2">
                        <Label htmlFor="email">Email</Label>
                        <Input id="email" type="email" placeholder="researcher@university.edu" value={email} onChange={(e) => setEmail(e.target.value)} required autoComplete="email" />
                      </div>
                      <div className="grid gap-2">
                        <div className="flex items-center"><Label htmlFor="password">Password</Label><a href="#" className="ml-auto text-xs underline-offset-2 hover:underline text-muted-foreground">Forgot password?</a></div>
                        <div className="relative">
                          <Input id="password" type={showPassword ? "text" : "password"} placeholder="••••••••" value={password} onChange={(e) => setPassword(e.target.value)} required autoComplete="current-password" className="pr-10" />
                          <button type="button" onClick={() => setShowPassword(!showPassword)} className="absolute right-3 top-1/2 -translate-y-1/2 text-muted-foreground hover:text-foreground transition-colors" tabIndex={-1}>
                            {showPassword ? <IconEyeOff className="size-4" /> : <IconEye className="size-4" />}
                          </button>
                        </div>
                      </div>
                      <Button type="submit" className="w-full" disabled={loading}>{loading ? "Signing in..." : "Sign in"}</Button>
                    </div>
                    <div className="relative text-center text-sm after:absolute after:inset-0 after:top-1/2 after:z-0 after:flex after:items-center after:border-t after:border-border"><span className="relative z-10 bg-card px-2 text-muted-foreground">or</span></div>
                    <div className="text-center text-sm text-muted-foreground">Kozi AI is in private beta. <a href="https://kozi-ai.com/contact" className="underline underline-offset-4 hover:text-primary font-medium" target="_blank" rel="noopener noreferrer">Request access</a></div>
                  </div>
                </form>
                <div className="relative hidden bg-gradient-to-br from-primary via-primary to-primary/80 md:flex md:flex-col md:items-center md:justify-center md:p-10">
                  <div className="absolute top-8 right-8 size-32 rounded-full bg-primary-foreground/5" />
                  <div className="absolute bottom-12 left-6 size-20 rounded-full bg-primary-foreground/5" />
                  <div className="relative z-10 flex flex-col items-center text-center text-primary-foreground max-w-xs">
                    <div className="flex size-16 items-center justify-center rounded-2xl bg-primary-foreground/10 backdrop-blur-sm mb-8"><IconFlask className="size-8 text-primary-foreground" /></div>
                    <p className="text-base font-medium leading-relaxed text-primary-foreground/90">&ldquo;From pathogen to candidate epitopes in minutes, not weeks.&rdquo;</p>
                    <div className="mt-6 flex items-center gap-2 text-sm text-primary-foreground/60"><div className="size-1.5 rounded-full bg-primary-foreground/40" /><span>Kozi AI</span></div>
                  </div>
                </div>
              </CardContent>
            </Card>
            <p className="text-center text-xs text-muted-foreground px-6">By continuing, you agree to our <a href="https://kozi-ai.com/terms" className="underline underline-offset-4 hover:text-primary" target="_blank" rel="noopener noreferrer">Terms</a> and <a href="https://kozi-ai.com/privacy" className="underline underline-offset-4 hover:text-primary" target="_blank" rel="noopener noreferrer">Privacy Policy</a>.</p>
          </div>
        </div>
      </div>
    </>
  );
}