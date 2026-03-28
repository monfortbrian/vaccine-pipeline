"use client"

import { usePathname, useRouter } from "next/navigation"
import { type Icon } from "@tabler/icons-react"
import { Badge } from "@/components/ui/badge"
import {
  SidebarGroup,
  SidebarGroupContent,
  SidebarMenu,
  SidebarMenuButton,
  SidebarMenuItem,
} from "@/components/ui/sidebar"

export function NavMain({
  items,
}: {
  items: {
    title: string
    url: string
    icon?: Icon
    soon?: boolean
  }[]
}) {
  const pathname = usePathname()
  const router = useRouter()

  return (
    <SidebarGroup>
      <SidebarGroupContent>
        <SidebarMenu>
          {items.map((item) => (
            <SidebarMenuItem key={item.title}>
              <SidebarMenuButton
                tooltip={item.title}
                isActive={pathname === item.url || pathname.startsWith(item.url + "/")}
                onClick={() => !item.soon && item.url !== "#" && router.push(item.url)}
                className={item.soon ? "opacity-60 cursor-default" : ""}
              >
                {item.icon && <item.icon />}
                <span>{item.title}</span>
                {item.soon && <Badge variant="secondary" className="ml-auto text-[9px] px-1.5 py-0">Soon</Badge>}
              </SidebarMenuButton>
            </SidebarMenuItem>
          ))}
        </SidebarMenu>
      </SidebarGroupContent>
    </SidebarGroup>
  )
}
