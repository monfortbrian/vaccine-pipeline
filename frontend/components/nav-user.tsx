"use client"

import { useRouter } from "next/navigation"
import {
  IconDotsVertical,
  IconLogout,
  IconUserCircle,
  IconSettings,
} from "@tabler/icons-react"
import { useAuth } from "@/components/auth-provider"
import { Avatar, AvatarFallback } from "@/components/ui/avatar"
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuGroup,
  DropdownMenuItem,
  DropdownMenuSeparator,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu"
import {
  SidebarMenu,
  SidebarMenuItem,
  useSidebar,
} from "@/components/ui/sidebar"

export function NavUser({
  user,
}: {
  user: {
    name: string
    email: string
    avatar: string
  }
}) {
  const { isMobile } = useSidebar()
  const { signOut } = useAuth()
  const router = useRouter()

  const handleLogout = async () => {
    await signOut()
    router.push("/login")
  }

  return (
    <SidebarMenu>
      <SidebarMenuItem>
        <DropdownMenu>
          <DropdownMenuTrigger>
            <div className="flex w-full items-center gap-2 rounded-lg p-2 text-sm hover:bg-sidebar-accent cursor-pointer">
              <Avatar className="h-8 w-8 rounded-full">
                <AvatarFallback className="rounded-full text-xs font-medium">
                  {user.name?.slice(0, 2)?.toUpperCase() || "KZ"}
                </AvatarFallback>
              </Avatar>
              <div className="grid flex-1 text-left text-sm leading-tight group-data-[collapsible=icon]:hidden">
                <span className="truncate font-medium">{user.name}</span>
                <span className="truncate text-xs text-muted-foreground">{user.email}</span>
              </div>
              <IconDotsVertical className="ml-auto size-4 group-data-[collapsible=icon]:hidden" />
            </div>
          </DropdownMenuTrigger>
          <DropdownMenuContent
            className="w-56 rounded-lg"
            side={isMobile ? "bottom" : "right"}
            align="end"
            sideOffset={4}
          >
            <div className="flex items-center gap-2 px-2 py-2 text-sm">
              <Avatar className="h-8 w-8 rounded-full">
                <AvatarFallback className="rounded-full text-xs font-medium">
                  {user.name?.slice(0, 2)?.toUpperCase() || "KZ"}
                </AvatarFallback>
              </Avatar>
              <div className="grid flex-1 text-left text-sm leading-tight">
                <span className="truncate font-medium">{user.name}</span>
                <span className="truncate text-xs text-muted-foreground">{user.email}</span>
              </div>
            </div>
            <DropdownMenuSeparator />
            <DropdownMenuGroup>
              <DropdownMenuItem onClick={() => router.push("/settings")}>
                <IconUserCircle className="mr-2 size-4" />Account
              </DropdownMenuItem>
              <DropdownMenuItem onClick={() => router.push("/settings")}>
                <IconSettings className="mr-2 size-4" />Settings
              </DropdownMenuItem>
            </DropdownMenuGroup>
            <DropdownMenuSeparator />
            <DropdownMenuItem onClick={handleLogout}>
              <IconLogout className="mr-2 size-4" />Log out
            </DropdownMenuItem>
          </DropdownMenuContent>
        </DropdownMenu>
      </SidebarMenuItem>
    </SidebarMenu>
  )
}
