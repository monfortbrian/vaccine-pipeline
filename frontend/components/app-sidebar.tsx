"use client"

import * as React from "react"
import {
  IconFlask,
  IconPlayerPlay,
  IconHistory,
  IconDna,
  IconChartBar,
  IconSettings,
  IconHelp,
  IconSearch,
  IconMolecule,
  IconFlask2,
  IconShieldCheck,
} from "@tabler/icons-react"

import { useAuth } from "@/components/auth-provider"
import { NavMain } from "@/components/nav-main"
import { NavSecondary } from "@/components/nav-secondary"
import { NavUser } from "@/components/nav-user"
import {
  Sidebar,
  SidebarContent,
  SidebarFooter,
  SidebarHeader,
  SidebarMenu,
  SidebarMenuButton,
  SidebarMenuItem,
} from "@/components/ui/sidebar"
import { Badge } from "@/components/ui/badge"

const navMain = [
  { title: "Playground", url: "/playground", icon: IconPlayerPlay },
  { title: "Run History", url: "/history", icon: IconHistory },
  { title: "Molecular Design", url: "#", icon: IconDna, soon: true },
  { title: "Compound Library", url: "#", icon: IconFlask2, soon: true },
  { title: "Validation", url: "#", icon: IconShieldCheck, soon: true },
]

const navSecondary = [
  { title: "Settings", url: "/settings", icon: IconSettings },
  { title: "Help", url: "https://kozi-ai.com", icon: IconHelp },
  { title: "Search", url: "#", icon: IconSearch },
]

export function AppSidebar({ ...props }: React.ComponentProps<typeof Sidebar>) {
  const { user } = useAuth()

  const userData = {
    name: user?.user_metadata?.name || user?.email?.split("@")[0] || "User",
    email: user?.email || "",
    avatar: "",
  }

  return (
    <Sidebar collapsible="icon" {...props}>
      <SidebarHeader>
        <SidebarMenu>
          <SidebarMenuItem>
            <SidebarMenuButton
              className="data-[slot=sidebar-menu-button]:!p-1.5"
              onClick={() => window.location.href = "/playground"}
            >
              <IconFlask className="!size-5" />
              <span className="text-base font-semibold">Kozi AI</span>
            </SidebarMenuButton>
          </SidebarMenuItem>
        </SidebarMenu>
      </SidebarHeader>
      <SidebarContent>
        <NavMain items={navMain} />
        <NavSecondary items={navSecondary} className="mt-auto" />
      </SidebarContent>
      <SidebarFooter>
        <NavUser user={userData} />
      </SidebarFooter>
    </Sidebar>
  )
}
