"use client";

import React, { createContext, useContext, useEffect, useState } from "react";
import { supabase } from "@/lib/supabase";
import type { User, Session } from "@supabase/supabase-js";

interface AuthCtx {
  user: User | null;
  session: Session | null;
  loading: boolean;
  signIn: (email: string, password: string) => Promise<void>;
  signOut: () => Promise<void>;
  updateDisplayName: (name: string) => Promise<void>;
}

const Ctx = createContext<AuthCtx | undefined>(undefined);

export function AuthProvider({ children }: { children: React.ReactNode }) {
  const [user, setUser] = useState<User | null>(null);
  const [session, setSession] = useState<Session | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    supabase.auth.getSession().then(({ data: { session: s } }) => {
      setSession(s);
      setUser(s?.user ?? null);
      setLoading(false);
    });
    const { data: { subscription } } = supabase.auth.onAuthStateChange((_, s) => {
      setSession(s);
      setUser(s?.user ?? null);
      setLoading(false);
    });
    return () => subscription.unsubscribe();
  }, []);

  const signIn = async (email: string, password: string) => {
    // First check if user is in beta_users with access
    const { data: betaUser } = await supabase
      .from("beta_users")
      .select("access_granted, name, organization, role")
      .eq("email", email)
      .single();

    // If not in beta_users OR not granted — show access denied
    if (!betaUser || !betaUser.access_granted) {
      throw new Error("ACCESS_DENIED");
    }

    // Now try actual login
    const { error } = await supabase.auth.signInWithPassword({ email, password });
    if (error) throw error;

    // Sync metadata
    if (betaUser.name) {
      await supabase.auth.updateUser({
        data: { name: betaUser.name, organization: betaUser.organization, role: betaUser.role }
      });
    }
  };

  const signOut = async () => {
    await supabase.auth.signOut();
    setUser(null);
    setSession(null);
  };

  const updateDisplayName = async (name: string) => {
    const { error } = await supabase.auth.updateUser({ data: { name } });
    if (error) throw error;
    if (user?.email) {
      await supabase.from("beta_users").update({ name }).eq("email", user.email);
    }
    const { data } = await supabase.auth.getUser();
    if (data.user) setUser(data.user);
  };

  return (
    <Ctx.Provider value={{ user, session, loading, signIn, signOut, updateDisplayName }}>
      {children}
    </Ctx.Provider>
  );
}

export function useAuth() {
  const c = useContext(Ctx);
  if (!c) throw new Error("useAuth must be within AuthProvider");
  return c;
}