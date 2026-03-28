"use client";

import { IconBrandX, IconBrandLinkedin, IconMail, IconFlask, IconBrandTiktok, IconBrandGithub } from "@tabler/icons-react";
import { Button } from "@/components/ui/button";
import { Card, CardContent } from "@/components/ui/card";

const CEO_EMAIL = "hugorshema@gmail.com";
const EMAIL_SUBJECT = encodeURIComponent("Beta Access Request - Kozi AI Platform");
const EMAIL_BODY = encodeURIComponent(`Dear Shema,

I recently tried to access the Kozi AI platform and my account is pending approval.

A bit about my work:
- Name: [Your Name]
- Institution: [Your Institution]
- Role: [Researcher / Bioinformatician / Student / Other]
- Research focus: [Brief description of your work]

I'm particularly interested in using Kozi for:
[Describe your use case - e.g., TB vaccine candidate screening, malaria epitope prediction, etc.]

I'd appreciate the opportunity to join the beta program. Happy to schedule a brief call to discuss how Kozi could support our research.

Best regards,
[Your Name]`);

export function AccessDeniedModal({ onClose }: { onClose: () => void }) {
  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/60 backdrop-blur-sm p-4">
      <Card className="w-full max-w-md animate-in fade-in zoom-in-95 duration-200">
        <CardContent className="flex flex-col items-center text-center p-8 gap-6">
          <div className="flex size-14 items-center justify-center rounded-2xl bg-primary/10">
            <IconFlask className="size-7 text-primary" />
          </div>

          <div className="space-y-2">
            <h2 className="text-xl font-bold">Request Access</h2>
            <p className="text-sm text-muted-foreground leading-relaxed max-w-sm">
              Oops, you’re seeing this for a reason. All good 🙂

            </p>
            <p className="text-sm text-muted-foreground leading-relaxed max-w-sm">
              First off, welcome to Kozi.
              We&apos;re rethinking how Labs should work.
            </p>
          </div>

          <div className="w-full rounded-lg bg-muted/50 p-4 text-left text-sm space-y-2">
            <p className="font-medium text-xs uppercase tracking-wider text-muted-foreground">What to do next?</p>
            <ul className="space-y-1.5 text-muted-foreground text-xs">
              <li>• You&apos;re trying to enter the playground, love that</li>
              <li>• No stress, we got you</li>
              <li>• Ping the CEO directly (yes, he reads them)</li>
              <li>• Get approved → jump right in 🚀</li>
            </ul>
          </div>

          <div className="w-full space-y-3">
            <p className="text-xs text-muted-foreground font-medium uppercase tracking-wider">
              Follow our journey
            </p>
            <div className="flex justify-center gap-3">
              <Button variant="outline" size="icon" className="rounded-full" onClick={() => window.open("https://x.com/kozi_ai", "_blank")}>
                <IconBrandX className="size-4" />
              </Button>
              <Button variant="outline" size="icon" className="rounded-full" onClick={() => window.open("https://linkedin.com/company/kozi-ai", "_blank")}>
                <IconBrandLinkedin className="size-4" />
              </Button>
              <Button variant="outline" size="icon" className="rounded-full" onClick={() => window.open("https://tiktok.com/@kozi_ai", "_blank")}>
                <IconBrandTiktok className="size-4" />
              </Button>
            </div>
          </div>

          <div className="w-full space-y-2">
            <Button className="w-full" onClick={() => window.open(`mailto:${CEO_EMAIL}?subject=${EMAIL_SUBJECT}&body=${EMAIL_BODY}`, "_blank")}>
              <IconMail className="size-4 mr-2" />
              Ping the CEO
            </Button>
            <Button variant="outline" className="w-full" onClick={() => window.open("https://kozi-ai.com/contact", "_blank")}>
              Request Access on Website
            </Button>
            <Button variant="ghost" className="w-full text-muted-foreground" onClick={onClose}>
              Back to login
            </Button>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}